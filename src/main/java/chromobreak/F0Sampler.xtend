package chromobreak

import bayonet.distributions.Random
import blang.mcmc.SampledVariable
import blang.mcmc.ConnectedFactor
import blang.core.Factor
import java.util.List
import blang.mcmc.MHSampler
import blang.mcmc.internals.Callback
import blang.types.internals.RealScalar

class F0Sampler extends MHSampler {
  
  @SampledVariable
  SingleCell model
  
  @ConnectedFactor
  List<Factor> factors
  
  override propose(Random random, Callback callback) {
    val f0 = model.f0 as RealScalar
    val oldValue = f0.doubleValue
    var ploidyJump = if (random.nextBoolean) Math::log(2) else Math::log(3)
    if (random.nextBoolean)
      ploidyJump *= -1
    val newValue = oldValue + ploidyJump
    f0.set(newValue)
    callback.proposalLogRatio = 0.0
    if (!callback.sampleAcceptance) 
      f0.set(oldValue)
  }
}