package humi.v8

import org.eclipse.xtend.lib.annotations.Data
import blang.core.IntDistribution
import blang.core.RealVar
import blang.inits.experiments.tabwriters.TidilySerializable
import blang.inits.experiments.tabwriters.TidySerializer.Context
import humi.HumiData
import blang.types.Index
import blang.types.Plated
import humi.CountFrequencies
import briefj.collections.Counter
import bayonet.math.NumericalUtils

class DistributionSummary implements TidilySerializable {
  var IntDistribution dist = null
  
  // The rate modelling the total count in the observed histogram, i.e. observed + unobserved
  var RealVar rate
  
  var int visCloneN_data
  var double visCloneTM_data
  
  val cutoff = 20
  
  def void init(IntDistribution dist, RealVar rate, CountFrequencies data) {
    this.dist = dist
    this.rate = rate
    visCloneN_data = data.nDataPoints
    visCloneTM_data = mean(truncatedNormalizedCounter(data))
  }
  
  def Counter<Integer> truncatedNormalizedCounter(CountFrequencies frequencies) {
    val result = new Counter
    for (c : 1 .. cutoff)
      if (frequencies.distinctCounts.contains(c))
        result.setCount(c, frequencies.frequency(c))
    result.normalize
    return result
  }
  
  def Counter<Integer> truncatedNormalizedCounter(IntDistribution dist) {
    val result = new Counter
    for (c : 1 .. cutoff)
      result.setCount(c, Math.exp(dist.logDensity(c)))
    result.normalize
    return result
  }
  
  // input: a normalized int valued dist, stored in a counter
  def double mean(Counter<Integer> counter) {
    NumericalUtils::checkIsClose(1.0, counter.totalCount) 
    var sum = 0.0
    for (c : counter.keySet)
      sum += c * counter.getCount(c)
    return sum
  }
  
  override serialize(Context context) {
    if (dist === null) throw new RuntimeException("Need to init first")
    val p0 = Math.exp(dist.logDensity(0))
    
    context.recurse((1.0 - p0) * rate.doubleValue, "description", "visCloneN_model")
    context.recurse(visCloneN_data,                "description", "visCloneN_data")
    
    context.recurse(mean(truncatedNormalizedCounter(dist)), "description", "visCloneTM_model")
    context.recurse(visCloneTM_data,                        "description", "visCloneTM_data")
  }
}