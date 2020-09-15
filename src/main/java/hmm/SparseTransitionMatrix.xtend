package hmm

import org.eclipse.xtend.lib.annotations.Data
import org.eclipse.xtend.lib.annotations.Delegate
import org.eclipse.xtend.lib.annotations.Accessors
import blang.types.internals.Delegator
import xlinear.SparseMatrix

import blang.types.TransitionMatrix

/** Matrix where each row is a SparseSimplex. */
@Data 
class SparseTransitionMatrix implements TransitionMatrix, SparseMatrix, Delegator<SparseMatrix> {
  @Accessors(PUBLIC_GETTER)
  @Delegate
  val SparseMatrix delegate 
  
  /** Get a view into a row. */
  override SparseSimplex row(int i) {
    return new SparseSimplex(delegate.row(i)) 
  }

}