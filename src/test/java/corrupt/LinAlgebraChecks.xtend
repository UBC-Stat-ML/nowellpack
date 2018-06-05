package corrupt

import org.junit.Test
import blang.validation.DiscreteMCTest
import java.util.ArrayList

import static extension corrupt.CorruptExtensionUtils.*
import static corrupt.CorruptStaticUtils.*

import static corrupt.EnumerationUtils.enumerate
import blang.runtime.SampledModel

class LinAlgebraChecks {
  
  val static nCells = 3
  val static nLoci = 2
  
  @Test
  def void test() {
    val list = enumerate(syntheticCells(nCells), syntheticLoci(nLoci))
    val equality = [SampledModel m | (m.model as Uniform).phylo.tree]
    val test = new DiscreteMCTest(list, equality)
    test.checkInvariance
    test.checkIrreducibility
  }
}