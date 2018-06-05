package corrupt

import org.junit.Test
import blang.validation.Instance
import blang.validation.DeterminismTest

import static extension briefj.BriefCollections.pick
import static extension corrupt.CorruptExtensionUtils.*
import static corrupt.CorruptStaticUtils.*

class EITest {
  
  val refCell = syntheticCells(1).pick
  val refLocus = syntheticLoci(1).pick
  
  val testFunction = [Uniform model |
    val tipIndic = model.phylo.tipIndicator(refCell, refLocus) 
    if (tipIndic.isIncluded) 1.0 else 0.0
  ]
  
  @Test
  def void testDerminism() {
    val phylo = new PerfectPhylo(syntheticCells(10), syntheticLoci(5))
    val model = new Uniform.Builder().setPhylo(phylo).build
    new DeterminismTest => [
      check(new Instance(model, testFunction))
    ]
  }
}