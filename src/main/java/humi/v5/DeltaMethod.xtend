package humi.v5

import binc.Command
import blang.inits.Arg
import blang.inits.DefaultValue
import blang.inits.experiments.Experiment
import blang.inits.parsing.Posix
import blang.types.Index
import blang.types.Plated
import briefj.BriefIO
import humi.CountFrequencies
import humi.HumiData
import humi.HumiStaticUtils
import java.util.ArrayList
import java.util.Collections

class DeltaMethod extends Experiment {
  @Arg HumiData data
  @Arg Plated<Integer> initialPopCounts
  
  @Arg                  @DefaultValue("0.95")
  public double winsorizedTailCutoff = 0.95
  
  /*
   * Value determined empirically to achieve approximate 95% coverage.
   * See ./nextflow run replicate-calibration.nf -resume --critical 5
   * got: 96% coverage between X0569_422234{1,2} and {2,3} and {1,3}
   * tried --critical 4 and got as low as 88% coverage. 
   * Theoretical value of 1.96 works on synthetic data (see bootstrap-calibration) 
   * but fails on replicate-calibration (as low as 60% coverage).
   */
  @Arg           @DefaultValue("5")
  public double criticalValue = 5  
  
  @Arg   @DefaultValue("Rscript")
  public String rCmd = "Rscript"
  
  override run() {
    compute(true)
    plot
  }
  
  static enum Columns { logRatio, logRatioIntervalRadius }
  
  def void compute(boolean withIntervals) {
    val controlInitialCounts = controlInitialCounts
    for (experiment : data.experiments.indices) {
      val writer = results.getTabularWriter("estimates").child(experiment.plate.name, experiment.key)
      val controlFinalCounts = controlFinalCounts(experiment)
      val controlFGRatio = controlFGRatio(experiment)
      for (gene : data.genes.indices)
        for (sgRNA : data.targets.indices(gene)) {
          val targetInitialCounts = initialPopCounts.get(sgRNA)
          val targetFinalCounts = finalCounts(sgRNA, experiment, false)
          val estimate = (targetFinalCounts / controlFinalCounts) / (targetInitialCounts / controlInitialCounts)
          val fgRatio = finalCounts(sgRNA, experiment, true) / (finalCounts(sgRNA, experiment,false) ** 2)
          val interval = criticalValue * Math::sqrt(controlFGRatio + fgRatio + 1.0/controlInitialCounts + 1.0/targetInitialCounts)
          (if (withIntervals)  
            writer.child(Columns::logRatioIntervalRadius, interval)
          else writer
          ).write(
            sgRNA.plate.name -> sgRNA.key,
            gene.plate.name -> gene.key,
            Columns::logRatio -> Math::log(estimate)
          )
        }
    }
    results.flushAll
  }
  
  def plot() {
    val plotResults = results.child("plots")
    val scriptFile = plotResults.getFileInResultFolder("script.r")
    BriefIO::write(scriptFile, '''
      require("ggplot2")
      data <- read.csv("«results.getFileInResultFolder("estimates.csv").absolutePath»")
      cols = rainbow(200, s=.6, v=.9)[sample(1:200,200)]
      p <- ggplot(data, aes(x = factor(«data.genes.name»), y = «Columns::logRatio», colour = factor(«data.targets.name»))) + 
        coord_flip() + 
        geom_errorbar(aes(ymin=«Columns::logRatio» - «Columns::logRatioIntervalRadius», ymax=«Columns::logRatio» + «Columns::logRatioIntervalRadius»)) +
        geom_point() + 
        facet_grid(. ~ «data.experiments.name») + 
        theme_bw() + 
        xlab("Gene") + 
        ylab("log(ratio)") + 
        ggtitle("Ratio of clone sizes relative to controls", subtitle = "Delta-method intervals") + 
        scale_colour_manual(values=cols) + 
        geom_hline(yintercept=0) + 
        theme(legend.position="none") 
      ggsave("«plotResults.getFileInResultFolder("intevals.pdf").absolutePath»", height = 10)
    ''')
    Command.call(Command.cmd(rCmd).appendArg(scriptFile.getAbsolutePath()))
  }
  
  def controlFGRatio(Index<String> experiment) {
    var num = 0.0
    var denom = 0.0
    for (gene : data.genes.indices) {
      if (data.isControl(gene)) {
        for (sgRNA : data.targets.indices(gene)) {
          num += finalCounts(sgRNA, experiment, true)
          denom += finalCounts(sgRNA, experiment, false) ** 2
        }
      }
    }
    return num / denom
  }
  
  def controlFinalCounts(Index<String> experiment) {
    var sum = 0.0
    for (gene : data.genes.indices) {
      if (data.isControl(gene)) {
        for (sgRNA : data.targets.indices(gene)) {
          sum += finalCounts(sgRNA, experiment, false)
        }
      }
    }
    return sum
  }
  
  def finalCounts(Index<Integer> sgRNA, Index<String> experiment, boolean squared) {
    val hist = data.histograms.get(sgRNA, experiment)
    val _testFunction = conditionalWinsorizedMean(hist, winsorizedTailCutoff)
    val testFunction = if (squared) [_testFunction.apply(it) ** 2] else _testFunction
    return sum(hist, testFunction)
  }
  
  def controlInitialCounts() {
    var sum = 0.0
    for (gene : data.genes.indices) {
      if (data.isControl(gene)) {
        for (sgRNA : data.targets.indices(gene)) {
          sum += initialPopCounts.get(sgRNA)
        }
      }
    }
    return sum
  }
  
  def static (Integer)=>Double conditionalWinsorizedMean(CountFrequencies hist, double p) {
    if (p < 0.5 || p > 1.0) throw new RuntimeException
    val distinctCounts = new ArrayList(hist.distinctCounts)
    Collections.sort(distinctCounts)
    var sum = 0.0
    var i = -1
    val normalization = hist.nDataPoints as double
    while (sum/normalization < p) {
      i++
      sum += hist.frequency(distinctCounts.get(i))
    }
    val cutOff = distinctCounts.get(i)
    return [count | (if (count < cutOff) count else cutOff).doubleValue]
  }
  
  def static sum(CountFrequencies hist, (Integer)=>Double f) {
    if (f.apply(0) != 0.0) throw new RuntimeException
    var result = 0.0
    for (count : hist.distinctCounts)
      result += f.apply(count) * hist.frequency(count) 
    return result
  }
  
  def static void main(String [] args) {
    Experiment::start(args, Posix::parse(args), HumiStaticUtils::parsingConfigs)
  }
}