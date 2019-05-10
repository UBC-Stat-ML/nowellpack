package humi

import java.util.Map
import java.util.LinkedHashMap
import blang.mcmc.Sampler
import bayonet.distributions.Random
import blang.mcmc.Samplers
import blang.mcmc.SampledVariable
import java.util.Collection


@Samplers(DummySampler) // Used for testing purpose (for EI tests)
class SimpleCountFrequencies implements CountFrequencies {
  public val Map<Integer, Integer> data = new LinkedHashMap // count -> frequency
  
  override distinctCounts() {
    return data.keySet
  }
  override frequency(int count) {
    return data.get(count).intValue
  }
  override nDataPoints() { 
    return data.values.stream.mapToInt[intValue].sum
  }
  
  def void add(int count) {
    val updated = data.getOrDefault(count, 0) + 1
    data.put(count, updated)
  }
  
  def void addAll(Collection<Integer> counts) {
    for (c : counts)
      add(c)
  }
  
  def void dropZeros() {
    data.remove(0)
  }
  
  static class DummySampler implements Sampler {
    @SampledVariable public SimpleCountFrequencies scf
    override execute(Random rand) {}
  }

  override toString() { toString(this) }
}