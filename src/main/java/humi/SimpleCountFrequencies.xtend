package humi

import java.util.Map
import java.util.LinkedHashMap
import blang.mcmc.Sampler
import bayonet.distributions.Random
import blang.mcmc.Samplers
import blang.mcmc.SampledVariable

// Used for testing purpose (for EI tests)
@Samplers(DummySampler)
class SimpleCountFrequencies implements CountFrequencies {
  public val Map<Integer, Integer> data = new LinkedHashMap
  
  override counts() {
    return data.keySet
  }
  override frequency(int count) {
    return data.get(count).intValue
  }
  
  static class DummySampler implements Sampler {
    @SampledVariable SimpleCountFrequencies scf
    override execute(Random rand) {}
  }
}