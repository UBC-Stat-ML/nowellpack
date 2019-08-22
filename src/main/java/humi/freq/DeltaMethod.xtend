package humi.freq

import blang.inits.Arg
import blang.inits.DefaultValue
import blang.inits.experiments.Experiment
import blang.inits.parsing.Posix
import blang.types.Index
import blang.types.Plated
import humi.CountFrequencies
import humi.HumiData
import humi.HumiStaticUtils
import java.util.ArrayList
import java.util.Collections

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
    HumiStaticUtils::plotIntervals(results, data, rCmd, true, "Delta-method intervals")
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