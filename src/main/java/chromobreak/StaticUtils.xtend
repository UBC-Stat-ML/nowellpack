package chromobreak

import xlinear.MatrixOperations
import static extension xlinear.MatrixExtensions.*
import blang.types.DenseSimplex

class StaticUtils {
  
  /**
   * Truncate zero and one and everything after max
   */
  def static DenseSimplex truncatedGeometric(double p, int max) {
    if (p.doubleValue <= 0 || p.doubleValue >= 1) blang.types.StaticUtils::invalidParameter
    val matrix = MatrixOperations::dense(max + 1)
    for (i : 2 .. max) 
      matrix.set(i, p * Math::pow(1.0 - p, i))
    matrix /= matrix.sum
    blang.types.StaticUtils::fixedSimplex(matrix)
  }
}