package corrupt

import bayonet.distributions.Random
import blang.core.LogScaleFactor
import blang.mcmc.ConnectedFactor
import blang.mcmc.SampledVariable
import blang.mcmc.Sampler
import blang.mcmc.internals.ExponentiatedFactor
import blang.mcmc.internals.SamplerBuilderContext

//import blang.distributions.Generators;
//import java.util.Random;

import blang.inits.Arg
import blang.inits.DefaultValue



class CorruptGibbsSampler implements Sampler {
  @Arg  @DefaultValue("1")
  public int overSample = 1
  
  @SampledVariable public CorruptPhylo phylo
  @ConnectedFactor public LogScaleFactor numericFactor
  
  private def double annealParameter() {
    val cast = numericFactor as ExponentiatedFactor
    val annealParam = cast.annealingParameter
    if (annealParam === null) return 1.0
    return annealParam.doubleValue
  }
  
  override execute(Random rand) {
    if (useTest)
      phylo.gibbsTest(rand, annealParameter)
    else {
	    	//val boolean resample = Random.nextBernoulli(new Random(10), skipProb) 
		// TODO: Monitor ESS and only re-sample if it has dropped below a certain threshold
		for (int p : 0 ..< overSample) {
			phylo.nextGibbs(rand, annealParameter)
		}
    }
  }
  
  override setup(SamplerBuilderContext context) {
    return true
  }
  
  public static var useTest = false // only for testing
}

