package humi.v5

import blang.inits.Arg
import humi.HumiData
import blang.types.Plated
import blang.inits.experiments.Experiment
import blang.inits.experiments.ParsingConfigs
import blang.io.internals.GlobalDataSourceStore
import blang.io.Parsers
import blang.inits.providers.CoreProviders
import blang.inits.Creators
import blang.inits.parsing.Posix
import blang.inits.parsing.Arguments
import blang.inits.Creator
import humi.CountFrequencies
import java.util.ArrayList
import java.util.Collections
import blang.types.Index
import blang.inits.DefaultValue
import briefj.BriefIO
import binc.Command

class DeltaMethod extends Experiment {
  @Arg HumiData data
  @Arg Plated<Integer> initialPopCounts
  
  @Arg @DefaultValue("0.95")
  public double winsorizedTailCutoff = 0.95
  
  @Arg   @DefaultValue("Rscript")
  public String rCmd = "Rscript"
  
  override run() {
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
          val interval = 1.96 * Math::sqrt(controlFGRatio + fgRatio + 1.0/controlInitialCounts + 1.0/targetInitialCounts)
          writer.write(
            sgRNA.plate.name -> sgRNA.key,
            gene.plate.name -> gene.key,
            "estimate" -> estimate,
            "logEstimate" -> Math::log(estimate),
            "leftInterval" -> Math::log(estimate) - interval,
            "rightInterval" -> Math::log(estimate) + interval,
            "width" -> interval
          )
        }
    }
    results.flushAll
    plot
  }
  
  def plot() {
    val plotResults = results.child("plots")
    val scriptFile = plotResults.getFileInResultFolder("script.r")
    BriefIO::write(scriptFile, '''
      require("ggplot2")
      data <- read.csv("«results.getFileInResultFolder("estimates.csv").absolutePath»")
      cols = rainbow(200, s=.6, v=.9)[sample(1:200,200)]
      p <- ggplot(data, aes(x = factor(gene), y = logEstimate, colour = factor(sgrna))) + 
        coord_flip() + 
        geom_errorbar(aes(ymin=leftInterval, ymax=rightInterval)) +
        geom_point() + 
        facet_grid(. ~ dataset) + 
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
    val Arguments parsedArgs = Posix.parse(args)
    val Creator creator = Creators::empty()
    creator.addFactories(CoreProviders)
    creator.addFactories(Parsers)
    val GlobalDataSourceStore globalDS = new GlobalDataSourceStore
    creator.addGlobal(GlobalDataSourceStore, globalDS)
    
    val ParsingConfigs parsingConfigs = new ParsingConfigs
    parsingConfigs.setCreator(creator) 
    Experiment::start(args, parsedArgs, parsingConfigs)
  }
}