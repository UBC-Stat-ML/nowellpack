package corrupt

import org.junit.Test

import static extension blang.engines.internals.MathCommonsAutoDiff.*

import xlinear.Matrix
import static  xlinear.MatrixOperations.*
import bayonet.distributions.Random
import corrupt.post.DeltaMethod.Transformation
import corrupt.post.DeltaMethod
import org.junit.Assert
import org.apache.commons.math3.stat.descriptive.SummaryStatistics

class TestDeltaMethod {
  
  val rand = new Random(1)
  
  val mu1 = 0.1
  val mu2 = -4.4
  val means = denseCopy(#[mu1, mu2])
  
  val rho = 0.17
  val var1 = 1.3
  val var2 = 5.1
  
  val cov = denseCopy(#[
    #[var1, rho * var1 * var2],
    #[rho * var1 * var2, var2]
  ])
  
  // using example where close form is easy, 
  // here https://www.stat.cmu.edu/~larry/=stat705/Lecture4.pdf
  
  val Transformation transform = [get(0) * get(1)]
  
  def Matrix generateData(int n) {
    val result = dense(n,2) 
    for (i : 0 ..< n) {
      val sample = sampleNormalByCovariance(rand, cov)
      sample += means
      result.row(i) += sample.transpose
    }
    return result
  }
  
  @Test
  def void testAsymptoticVariance() {
    val delta = new DeltaMethod(generateData(10000000), transform) 
    val numerical = delta.asymptoticVariance
    val analytical = mu2 * mu2 * cov.get(0,0) + 2.0 * mu1 * mu2 * cov.get(0,1) + mu1 * mu1 * cov.get(1,1)
    Assert.assertEquals(numerical, analytical, 0.1)
  }
  
  @Test
  def void testCoverage() {
    val coverage = 0.7
    val truth = mu1 * mu2
    val stats = new SummaryStatistics
    for (rep : 0 ..< 10000) {
      val delta = new DeltaMethod(generateData(1000), transform)
      val ciRadius = delta.confidenceIntervalRadius(coverage)
      val estimate = delta.estimate
      val means = delta.means
      Assert.assertEquals(estimate, means.get(0) * means.get(1), 1e-20)
      val contained = estimate - ciRadius <= truth && truth <= estimate + ciRadius
      stats.addValue(if (contained) 1.0 else 0.0)
    }
    Assert.assertEquals(coverage, stats.mean, 0.2)
  }
}