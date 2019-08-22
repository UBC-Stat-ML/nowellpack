package humi.freq

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
import blang.inits.experiments.ExperimentResults

class DeltaMethod extends Experiment {
  @Arg HumiData data
  @Arg Plated<Integer> initialPopCounts
  
  @Arg                  @DefaultValue("0.99")
  public double winsorizedTailCutoff = 0.99

  @Arg           @DefaultValue("1.96")
  public double criticalValue = 1.96  
  
  @Arg   @DefaultValue("Rscript")
  public String rCmd = "Rscript"
  
  override run() {
    compute
    plot(results, data, rCmd, true, "Delta-method intervals")
  }
  
  static enum Columns { logRatio, logRatioLeftBound,  logRatioRightBound}
  
  def void compute() {
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
          val logEstimate = Math::log(estimate)
          val fgRatio = finalCounts(sgRNA, experiment, true) / (finalCounts(sgRNA, experiment,false) ** 2)
          val interval = criticalValue * Math::sqrt(controlFGRatio + fgRatio + 1.0/controlInitialCounts + 1.0/targetInitialCounts)
          writer.write(
            Columns::logRatioLeftBound -> logEstimate-interval,    
            Columns::logRatioRightBound -> logEstimate+interval,
            sgRNA.plate.name -> sgRNA.key,
            gene.plate.name -> gene.key,
            Columns::logRatio -> logEstimate
          )
        }
    }
    results.flushAll
  }
  
  def static plot(ExperimentResults results, HumiData data, String rCmd, boolean expSpecific, String caption) {
    val plotResults = results.child("plots")
    val scriptFile = plotResults.getFileInResultFolder("script.r")
    BriefIO::write(scriptFile, '''
      require("ggplot2")
      data <- read.csv("«results.getFileInResultFolder("estimates.csv").absolutePath»")
      cols = rainbow(200, s=.6, v=.9)[sample(1:200,200)]
      p <- ggplot(data, aes(x = factor(«data.genes.name»), y = «Columns::logRatio», colour = factor(«data.targets.name»))) + 
        coord_flip() + 
        geom_errorbar(aes(ymin=«Columns::logRatioLeftBound», ymax=«Columns::logRatioRightBound»)) +
        geom_point() + 
        «IF expSpecific»facet_grid(. ~ «data.experiments.name») + «ENDIF»
        theme_bw() + 
        xlab("Gene") + 
        ylab("log(ratio)") + 
        ggtitle("Ratio of clone sizes relative to controls", subtitle = "«caption»") + 
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