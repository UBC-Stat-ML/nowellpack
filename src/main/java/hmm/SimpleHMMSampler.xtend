package hmm

import blang.mcmc.Sampler
import bayonet.distributions.Random
import blang.mcmc.internals.SamplerBuilderContext
import blang.mcmc.ConnectedFactor
import java.util.List
import blang.core.LogScaleFactor
import blang.mcmc.SampledVariable
import blang.core.RealVar
import blang.core.WritableIntVar
import blang.core.Factor

class SimpleHMMSampler implements Sampler {
  
  @ConnectedFactor
  List<Factor> factors
  
  @SampledVariable
  SimpleHMM hmm
  
  RealVar anneal
  
  override execute(Random rand) {
    if (anneal.doubleValue != 1.0) throw new RuntimeException
    val dynamics = hmm.dynamics
    val emissions = hmm.emissions
    val initial = hmm.initial
    val obs = hmm.observations
    val dp = new HMM() {
      override transitionProbabilities(int t) {
        return dynamics
      }
      override initialProbabilities() {
        return initial
      }
      override length() {
        hmm.observations.size
      }
      override observationLogDensity(int t, int state) {
        Math::log(emissions.get(state, obs.get(t).intValue))
      }
    }
    val list = HMMComputations::sample(rand, dp)
    for (i : 0 ..< obs.size) {
      val cur = hmm.observations.get(i) as WritableIntVar
      cur.set(list.get(i))
    }
  }
  
  override boolean setup(SamplerBuilderContext context) {
    this.anneal = context.annealingParameter
    return true
  }  
}