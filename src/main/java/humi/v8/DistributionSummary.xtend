package humi.v8

import blang.core.IntDistribution
import blang.types.Index
import humi.CountFrequencies
import briefj.collections.Counter
import bayonet.math.NumericalUtils
import blang.types.internals.SimplePlate
import blang.types.Plate

class DistributionSummary {
  
  static String dataString = "data"
  static String modelString = "model"
  public static Plate<String> description = new SimplePlate("description", #{dataString, modelString})
  public static Index<String> dataIndex = new Index(description, dataString)
  public static Index<String> modelIndex = new Index(description, modelString)
  
  static val cutoff = 20
  
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