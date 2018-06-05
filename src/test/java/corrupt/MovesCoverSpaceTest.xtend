package corrupt

import org.junit.Test

import static java.lang.Math.exp
import bayonet.distributions.Random
import static corrupt.EnumerationUtils.enumerate

import static org.junit.Assert.assertEquals
import bayonet.math.NumericalUtils

import static extension corrupt.CorruptExtensionUtils.*
import static corrupt.CorruptStaticUtils.*
import java.util.ArrayList

class MovesCoverSpaceTest {
  
  val static nCells = 2
  val static nLoci = 3
  
  def void sampleNonUniform(Random random) {
    new PerfectPhylo(syntheticCells(nCells), syntheticLoci(nLoci)) => [
      for (split : splits.values) 
        split.tree.collapseEdge(split.locus)
      for (split : splits.values) {
        val parent = random.uniformElement(tree.lociAndRoot)
        val movedChildren = random.uniformSubset(new ArrayList(tree.children(parent)))
        split.tree.addEdge(parent, split.locus, movedChildren) 
      }
    ]
  }
  
  @Test
  def void test() {
    assertEquals(
      println(exp(logNPerfectPhylo(nCells, nLoci))), 
      enumerate(syntheticCells(nCells), syntheticLoci(nLoci)).size,  
      NumericalUtils.THRESHOLD
    )
  }
}