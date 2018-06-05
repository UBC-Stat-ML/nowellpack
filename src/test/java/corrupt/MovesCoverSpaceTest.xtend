package corrupt

import org.junit.Test

import static java.lang.Math.exp
import bayonet.distributions.Random

import static org.junit.Assert.assertEquals
import bayonet.math.NumericalUtils

import static extension corrupt.CorruptExtensionUtils.*
import static corrupt.CorruptStaticUtils.*
import java.util.ArrayList
import static corrupt.EnumerationUtils.enumerateUniformModels

class MovesCoverSpaceTest {
  
  val static nCells = 2
  val static nLoci = 3
  
  def void sampleNonUniform(Random random) {
    new PerfectPhylo(syntheticCells(nCells), syntheticLoci(nLoci)) => [
      for (split : splits.values) 
        split.tree.collapseEdge(split.locus)
      for (split : splits.values) {
        val parent = random.uniformElement(tree.lociAndRoot)
        val movedChildren = random.uniformSubset((tree.children(parent)))
        split.tree.addEdge(parent, split.locus, movedChildren) 
      }
    ]
  }
  
  @Test
  def void test() {
    assertEquals(
      println(exp(logNPerfectPhylo(nCells, nLoci))), 
      enumerateUniformModels(nCells, nLoci).size,  
      NumericalUtils.THRESHOLD
    )
  }
}