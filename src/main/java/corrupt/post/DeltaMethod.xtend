package corrupt.post

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure
import static extension blang.engines.internals.MathCommonsAutoDiff.*


import java.util.List
import org.apache.commons.math3.stat.correlation.Covariance
import xlinear.Matrix
import static extension xlinear.MatrixExtensions.*
import static  xlinear.MatrixOperations.*
import xlinear.internals.CommonsDenseMatrix
import org.apache.commons.math3.stat.descriptive.moment.Mean
import org.eclipse.xtend.lib.annotations.Data
import java.util.ArrayList
import org.apache.commons.math3.distribution.NormalDistribution

@Data
class DeltaMethod {
   
  val Matrix data 
  val Transformation transformation
  
 /**
   * Build a 95% confidence interval
   */
  def double confidenceIntervalRadius() {
    confidenceIntervalRadius(0.95)
  }
  
  /**
   * Construct the radius of an interval I such that asymptotically in n,
   * P(Z \in I) = coverage
   */
  def double confidenceIntervalRadius(double coverage) {
    if (coverage < 0.5 || coverage > 1.0)  // putting 0.5 instead of 1.0 to catch misinterpretations of input argument which is not alpha
      throw new RuntimeException
    val alpha = 1.0 - coverage
    val stdNormal = new NormalDistribution(0.0, 1.0)
    val z = stdNormal.inverseCumulativeProbability(1.0 - alpha/2.0)
    return z * Math.sqrt(asymptoticVariance / n)
  }
  
  /**
   * n = number of observations
   * p = number of parameters
   * Input: nxp matrix
   * Output: pxp covar matrix
   */
  def static Matrix covar(Matrix data) {
    return new CommonsDenseMatrix(new Covariance(data.toCommonsMatrix).covarianceMatrix)
  }
  
  def int n() { return data.nRows }
  def int p() { return data.nCols }
  
  def Matrix means() {
    val result = dense(p)
    for (i : 0 ..< p) {
      val m = new Mean()
      result.set(i, m.evaluate(data.col(i).vectorToArray))
    }
    return result
  }
  
  def double estimate() {
    derivativeStructure.value
  }
  
  def DerivativeStructure derivativeStructure() {
    val means = means()
    val args = new ArrayList
    for (i : 0 ..< p) {
      args.add(new DerivativeStructure(means.nEntries, 1, i, means.get(i)))
    }
    return transformation.apply(args)
  }
  
  def double asymptoticVariance() {
    val gradient = derivativeStructure.gradient
    val result = gradient.transpose * covar(data) * gradient
    return result.doubleValue
  }
  
  @FunctionalInterface
  static interface Transformation {
    def DerivativeStructure apply(List<DerivativeStructure> arguments)
  }
}