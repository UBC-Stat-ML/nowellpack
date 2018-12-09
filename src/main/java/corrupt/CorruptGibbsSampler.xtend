package corrupt

import bayonet.distributions.Random
import blang.core.LogScaleFactor
import blang.mcmc.ConnectedFactor
import blang.mcmc.SampledVariable
import blang.mcmc.Sampler
import blang.mcmc.internals.ExponentiatedFactor
import blang.mcmc.internals.SamplerBuilderContext
import java.util.Collections
import java.util.ArrayList

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
    if (useTest) 
      executeTest(rand)
    else {
      if (SamplerOptions::instance.useCellReallocationMove)
        phylo.sample(rand, annealParameter)
      else
        phylo.sampleWithoutCellReallocation(rand, annealParameter)
    }
  }
  
  def executeTest(Random rand) {
    if (rand.nextBernoulli(0.5))
      phylo.sample(rand, annealParameter, Collections::emptyList, new ArrayList => [add(phylo.cells.get(rand.nextInt(phylo.cells.size)))])
    else
      phylo.sample(rand, annealParameter, new ArrayList => [add(phylo.loci .get(rand.nextInt(phylo.loci .size)))], Collections::emptyList)
  }
  
  override setup(SamplerBuilderContext context) {
    return true
  }
  
  public static var useTest = false
}