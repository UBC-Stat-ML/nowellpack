package corrupt

import org.junit.Test
import static corrupt.CorruptUtils.syntheticLoci
import static corrupt.CorruptUtils.syntheticCells

import static extension corrupt.CorruptUtils.uniformElement
import static extension corrupt.CorruptUtils.uniformSubset
import static extension corrupt.CorruptUtils.lociAndRoot
import static corrupt.CorruptUtils.logNPerfectPhylo
import static java.lang.Math.exp
import bayonet.distributions.Random
import bayonet.distributions.ExhaustiveDebugRandom


import static org.junit.Assert.assertEquals
import bayonet.math.NumericalUtils

class MovesCoverSpaceTest {
  
  val static nCells = 3
  val static nLoci = 2
  
  def void sampleNonUniform(Random random) {
    new PerfectPhylo(syntheticCells(nCells), syntheticLoci(nLoci)) => [
      for (split : splits.values) 
        split.tree.collapseEdge(split.locus)
      for (split : splits.values) {
        val parent = random.uniformElement(tree.lociAndRoot)
        val movedChildren = random.uniformSubset(tree.children(parent))
        split.tree.addEdge(parent, split.locus, movedChildren) 
      }
    ]
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