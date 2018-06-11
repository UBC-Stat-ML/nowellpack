package corrupt

import bayonet.distributions.Random
import blang.core.LogScaleFactor
import blang.core.RealVar
import blang.mcmc.ConnectedFactor
import blang.mcmc.SampledVariable
import blang.mcmc.Sampler
import blang.mcmc.internals.ExponentiatedFactor
import blang.mcmc.internals.SamplerBuilderContext

class CorruptGibbsSampler implements Sampler {
  @SampledVariable public CorruptPhylo phylo
  @ConnectedFactor public LogScaleFactor numericFactor 
  
  private def double annealParameter() {
    val cast = numericFactor as ExponentiatedFactor
    val annealParam = cast.annealingParameter
    if (annealParam === null) return 1.0
    return annealParam.doubleValue
  }
  
  override execute(Random rand) {
    phylo.gibbSample(rand, annealParameter)
  }
  
  override setup(SamplerBuilderContext context) {
    
    return true
  }
}
