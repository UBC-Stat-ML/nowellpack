package corrupt

import org.junit.Test
import bayonet.distributions.Random
import corrupt.post.CLMatrixUtils
import corrupt.post.ReadOnlyCLMatrix
import org.junit.Assert

class TestGreedy {
 
  @Test
  def void test() {
    val nCells = 100
    val nLoci = 50
    val rand = new Random(1)
    for (i : 0 .. 10) {
      val phylo = PerfectPhylo::generateUniform(nCells, nLoci, rand)
      
      val matrix = ReadOnlyCLMatrix.readOnly(CLMatrixUtils::fromPhylo(phylo))  
      val greedy = new Greedy => [
        tipInclusionProbabilities = matrix
      ]
      val inferred = greedy.infer
      Assert::assertEquals(0, CLMatrixUtils::distance(phylo, inferred.reconstruction), 0)
    }
  }
}