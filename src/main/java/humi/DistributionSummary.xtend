package humi

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
  static val winsorizationP = 0.99
  
  def static registerMonitors(
    Plated<Monitor> visibleCloneNumbers, 
    Plated<Monitor> truncatedMeans, 
    Plated<Monitor> winsorizedMeans,
    Plated<Monitor> conditionalWinsorizedMeans,
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
    
    // winsorized mean
    winsorizedMeans.get(target).init(new RealVar() {
      override doubleValue() {
        return winsorizedMean(dist.get, winsorizationP)
      }
    })
    
    // conditional version 
    conditionalWinsorizedMeans.get(target).init(new RealVar() {
      override doubleValue() {
        return conditionalWinsorizedMean(dist.get, winsorizationP)
      }
    })
  }
  
  def static double conditionalWinsorizedMean(IntDistribution distribution, double p) {
    if (p < 0.5 || p > 1.0) throw new RuntimeException
    var sum = 0.0
    var x = 0
    val normalization = 1.0 - Math::exp(distribution.logDensity(0))
    while (sum/normalization < p) {
      x++  
      sum += Math::exp(distribution.logDensity(x))
    }
    val cutOff = x
    var result = 0.0
    var mass = 0.0
    for (y : 0 ..< cutOff) {
      val currentMass = Math::exp(distribution.logDensity(y))
      mass += currentMass
      result += y * currentMass
    }
    result += cutOff * (1.0 - mass)
    return result
  }
  
  def static double winsorizedMean(IntDistribution distribution, double p) {
    return mean(truncate(distribution, p))
  }
  
  def private static Counter<Integer> truncate(IntDistribution distribution, double p) {
    val result = new Counter
    var mass = 0.0
    var c = 0
    while (true) {
      val cur =  Math::exp(distribution.logDensity(c))
      
      if (mass + cur < p) {
        result.setCount(c, cur)
        mass += cur
      } else {
        result.setCount(c, 1.0 - mass)
        NumericalUtils::checkIsClose(result.totalCount, 1.0)
        return result
      }
      c++
    }
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