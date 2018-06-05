package corrupt

import org.junit.Test
import blang.validation.Instance
import blang.validation.DeterminismTest

import static extension briefj.BriefCollections.pick
import static extension corrupt.CorruptExtensionUtils.*
import static corrupt.CorruptStaticUtils.*
import blang.validation.ExactInvarianceTest

class EITest {
  
  val refCell = syntheticCells(1).pick
  val refLocus = syntheticLoci(1).pick
  
  val phylo = new PerfectPhylo(syntheticCells(10), syntheticLoci(5))
  val model = new Uniform.Builder().setPhylo(phylo).build
  
  val uniformInstance = new Instance(
    model, 
    [Uniform model |
      val tipIndic = model.phylo.tipIndicator(refCell, refLocus) 
      if (tipIndic.isIncluded) 1.0 else 0.0
    ]
  )
  
  val syntheticInstance = new Instance(
    EnumerationUtils::syntheticModel(8, 9), 
    [Synthetic model |
      val tipIndic = model.phylo.tipIndicator(refCell, refLocus) 
      if (tipIndic.isIncluded) 1.0 else 0.0
    ],
    [Synthetic model |
      model.observations.entries.iterator.next.value.doubleValue
    ] 
  )
  
  val instances = #[uniformInstance, syntheticInstance]
  
  @Test
  def void testDerminism() {
    
    new DeterminismTest => [
      for (i : instances)
        check(i)
    ]
  }
  
  @Test
  def void eit() {
    val test = new ExactInvarianceTest => [
      nPosteriorSamplesPerIndep = 500
      for (i : instances)
        add(i)
    ]
    test.check
  }
}