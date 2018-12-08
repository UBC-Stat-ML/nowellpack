package corrupt

import org.junit.Test
import blang.validation.DiscreteMCTest

import static extension corrupt.CorruptExtensionUtils.*
import static corrupt.CorruptStaticUtils.*

import static java.lang.Math.*

import blang.runtime.SampledModel
import static corrupt.EnumerationUtils.*

class LinAlgebraChecks {
  
  val static nCells = 3
  val static nLoci = 2
  
  @Test
  def void testMove() {
    CorruptGibbsSampler::useTest = true
    for (anneal : #[0.0, 0.42, 1.0]) {
      val list = enumerateSyntheticModels(nCells, nLoci, anneal) 
      val equality = [SampledModel m | (m.model as FixedMatrixModel).phylo.getReconstruction.tree] 
      val test = new DiscreteMCTest(list, equality, false)
      test.verbose = true
      test.checkInvariance
      test.checkIrreducibility
      test.checkStateSpaceSize(round(exp(logNPerfectPhylo(nCells, nLoci))) as int)
    }
  }
}