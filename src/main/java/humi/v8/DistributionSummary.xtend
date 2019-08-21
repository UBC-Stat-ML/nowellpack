package humi.v8

import blang.core.IntDistribution
import blang.types.Index
import humi.CountFrequencies
import briefj.collections.Counter
import bayonet.math.NumericalUtils
import blang.types.internals.SimplePlate
import blang.types.Plate
import humi.Monitor
import blang.types.Plated
import java.util.function.Supplier
import blang.core.IntVar
import blang.core.RealVar
import humi.HumiData
import blang.core.RealConstant

class DistributionSummary {
  
  static String dataString = "data"
  static String modelString = "model"
  public static Plate<String> description = new SimplePlate("description", #{dataString, modelString})
  public static Index<String> dataIndex = new Index(description, dataString)
  public static Index<String> modelIndex = new Index(description, modelString)
  
  static val cutoff = 20
  
  def static registerMonitors(
    Plated<Monitor> visibleCloneNumbers, 
    Plated<Monitor> truncatedMeans, 
    Index<Integer> target,
    Supplier<IntDistribution> dist,
    Plated<IntVar> initialPopCounts,
    RealVar lambda,
    HumiData data,
    Index<String> experiment
  ) {
    val initialPopCount = initialPopCounts.get(target)
    
    // register visible clone number monitors
    visibleCloneNumbers.get(target, DistributionSummary::modelIndex).init(new RealVar() {
      override doubleValue() {
        val p0 = Math.exp(dist.get.logDensity(0))
        return (1.0 - p0) * initialPopCount.intValue * lambda.doubleValue
      }
    })
    val observedHist = data.histograms.get(target, experiment)
    val dataVisibleCloneNumbers = observedHist.nDataPoints
    visibleCloneNumbers.get(target, DistributionSummary::dataIndex).init(new RealConstant(dataVisibleCloneNumbers))
    
    // register truncated mean monitors
    truncatedMeans.get(target, DistributionSummary::modelIndex).init(new RealVar() {
      override doubleValue() {
        DistributionSummary::mean(DistributionSummary::truncatedNormalizedCounter(dist.get))
      }
    })
    val dataCloneTM = DistributionSummary::mean(DistributionSummary::truncatedNormalizedCounter(observedHist))
    truncatedMeans.get(target, DistributionSummary::dataIndex).init(new RealConstant(dataCloneTM))
  }
  
  def static Counter<Integer> truncatedNormalizedCounter(CountFrequencies frequencies) {
    val result = new Counter
    for (c : 1 .. cutoff)
      if (frequencies.distinctCounts.contains(c))
        result.setCount(c, frequencies.frequency(c))
    result.normalize
    return result
  }
  
  def static Counter<Integer> truncatedNormalizedCounter(IntDistribution dist) {
    val result = new Counter
    for (c : 1 .. cutoff)
      result.setCount(c, Math.exp(dist.logDensity(c)))
    result.normalize
    return result
  }
  
  // input: a normalized int valued dist, stored in a counter
  def static double mean(Counter<Integer> counter) {
    NumericalUtils::checkIsClose(1.0, counter.totalCount) 
    var sum = 0.0
    for (c : counter.keySet)
      sum += c * counter.getCount(c)
    return sum
  }
}