package corrupt

import org.junit.Test
import static corrupt.CorruptUtils.syntheticLoci
import blang.validation.Instance
import blang.validation.DeterminismTest

class EITest {
  
  val testFunction = [Uniform model |
    val tipIndic = model.phylo.tipIndicator(new Cell(0) , syntheticLoci(1).iterator.next)
    if (tipIndic.isIncluded) 1.0 else 0.0
  ]
  
  @Test
  def void testDerminism() {
    val phylo = new PerfectPhylo(1, syntheticLoci(1))
    val model = new Uniform.Builder().setPhylo(phylo).build
    new DeterminismTest => [
      check(new Instance(model, testFunction))
    ]
  }
}