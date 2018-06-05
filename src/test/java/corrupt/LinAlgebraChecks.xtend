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
  def void test() {
    for (useData : #[false, true]) {
      val list = 
        if (useData)
          enumerateSyntheticModels(nCells, nLoci) 
        else
          enumerateUniformModels(nCells, nLoci)
      val equality = 
        if (useData) 
          [SampledModel m | (m.model as Synthetic).phylo.tree]
        else
          [SampledModel m | (m.model as Uniform).phylo.tree]
      val test = new DiscreteMCTest(list, equality, false)
      test.verbose = true
      test.checkInvariance
      test.checkIrreducibility
      test.checkStateSpaceSize(round(exp(logNPerfectPhylo(nCells, nLoci))) as int)
    }
  }
}