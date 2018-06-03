package corrupt

import org.junit.Test
import static corrupt.CorruptUtils.syntheticLoci

import static extension corrupt.CorruptUtils.uniformElement
import static extension corrupt.CorruptUtils.uniformSubset
import static extension corrupt.CorruptUtils.lociAndRoot
import static corrupt.CorruptUtils.logNPerfectPhylo
import static java.lang.Math.exp
import bayonet.distributions.Random
import bayonet.distributions.ExhaustiveDebugRandom


import static org.junit.Assert.assertEquals
import bayonet.math.NumericalUtils

class MovesCoverSpaceTest extends PerfectPhylo {
  
  val static nCells = 3
  val static nLoci = 2
  
  new() {
    super(nCells, syntheticLoci(nLoci))
  }
  
  def void sampleNonUniform(Random random) {
    for (split : splits.values) 
      split.remove
    for (split : splits.values) {
      val parent = random.uniformElement(tree.lociAndRoot)
      val movedChildren = random.uniformSubset(tree.children(parent))
      split.reAttach(parent, movedChildren) 
    }
  }
  
  @Test
  def void test() {
    val exhaustive = new ExhaustiveDebugRandom
    var count = 0
    while (exhaustive.hasNext) {
      sampleNonUniform(exhaustive)
      count++
    }
    assertEquals(exp(logNPerfectPhylo(nCells, nLoci)), count, NumericalUtils.THRESHOLD)
  }
}