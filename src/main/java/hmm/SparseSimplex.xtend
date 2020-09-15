package hmm

import org.eclipse.xtend.lib.annotations.Delegate
import blang.mcmc.Samplers
import blang.mcmc.SimplexSampler
import org.eclipse.xtend.lib.annotations.Accessors
import blang.types.internals.Delegator
import xlinear.SparseMatrix
import bayonet.math.NumericalUtils

import static extension xlinear.MatrixExtensions.*
import blang.types.Simplex

/** Vector of entries summing to one.
 * 
 * We do not enforce positive constraints here to facilitate undoing sampling moves.
 */
@Samplers(SimplexSampler)
class SparseSimplex implements Simplex, SparseMatrix, Delegator<SparseMatrix> {
  @Accessors(PUBLIC_GETTER)
  @Delegate SparseMatrix delegate
  
  new (SparseMatrix matrix) {
    NumericalUtils::checkIsClose(matrix.sum, 1.0)
    this.delegate = matrix
  }
  
  /** Set a pair of entries, checking their sum is the same before and after */
  def void setPair(int index1, double value1, int index2, double value2) {
    val double old = get(index1) + get(index2)
    NumericalUtils::checkIsClose(old, value1 + value2)
    delegate.set(index1, value1)
    delegate.set(index2, value2)
  }
  
  override void set(int i, int j, double value) {
    throw new RuntimeException("Use setPair instead");
  }
  
  override void set(int i, double value) {
    throw new RuntimeException("Use setPair instead");
  }
  
}