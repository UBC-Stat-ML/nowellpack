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
import java.util.List

class CorruptGibbsSampler implements Sampler {
  @SampledVariable public CorruptPhylo phylo
  @ConnectedFactor public LogScaleFactor numericFactor
  
  private def double annealParameter() {
    val cast = numericFactor as ExponentiatedFactor
    val annealParam = cast.annealingParameter
    if (annealParam === null) return 1.0
    return annealParam.doubleValue
  }
  
  var List<Locus> loci = null
  var List<Cell> cells = null
  var int i = 0
  
  def static <T> T get(List<T> items, int i, Random rand) {
    if (i === 0) Collections.shuffle(items, rand)
    return items.get(i % items.size)
  }
  
  override execute(Random rand) {
    if (useTest) 
      executeTest(rand)
    else {
      if (SamplerOptions::instance.useMiniMoves) {
        if (loci === null) loci = new ArrayList(phylo.loci)
        if (cells === null) cells = new ArrayList(phylo.cells)
        val sampledLoci = #[get(loci, i, rand)]
        val sampledCells = if (SamplerOptions::instance.useCellReallocationMove) #[get(cells, i, rand)] else #[]
        phylo.sample(rand, annealParameter, sampledLoci, sampledCells)
        i++
      } else {
        if (SamplerOptions::instance.useCellReallocationMove)
          phylo.sample(rand, annealParameter)
        else
          phylo.sampleWithoutCellReallocation(rand, annealParameter)
      }
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