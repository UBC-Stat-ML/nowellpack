package corrupt.post

import corrupt.Cell
import corrupt.Locus
import java.util.Set
import blang.inits.Implementations

/**
 * The primary usage of this interface is to represent, for 
 * each locus l and cell c, the likelihood P(Y_{cl} = y | X_{cl} = x), 
 * where Y_{cl} is the observation and X_{cl} the latent variable 
 * indicating whether cell c has trait l.
 * 
 * For the purpose of efficient recursion as well as a storage, we 
 * store only one number for each cell and trait. This is done as follows.
 * 
 * Assume a uniform pseudo-prior giving uniform iid probability to each X_{cl}. 
 * This pseudo-prior is purely an algorithmic intermediate quantity. 
 * Under this pseudo model, at each locus and cell, we obtain:
 * 
 * P(Y_{cl} = y | X_{cl} = x) = Z_{cl} P(X_{cl} = x | Y_{cl} = y)
 * Z_{cl} = P(Y_{cl} = y | X_{cl} = 0) + P(Y_{cl} = y | X_{cl} = 1)
 * 
 * The quantity we need in likelihood calculations is P(Y_{cl} = y | X_{cl} = x). 
 * 
 * Instead, we provide P(X_{cl} = 1 | Y_{cl} = y), from which algorithms can instantiate 
 * 
 * P(X_{cl} = 0 | Y_{cl} = y) = 1 - P(X_{cl} = 1 | Y_{cl} = y)
 * 
 * To get a valid likelihood in the context of parameter estimation, we also provide the log of the product 
 * log(Z) = log(\prod_c \prod_l Z_{cl})
 * 
 * 
 * More generally, other matrices indexed by cell and locus are also 
 * represented for convenience and simplicity using this data structure. In such case the normalization 
 * is set to one (i.e. log set to zero).
 */
@Implementations(ReadOnlyCLMatrix, NoisyBinaryCLMatrix)
interface CellLocusMatrix {
  /**
   * In the context of constructing a likelihood, P(X_{cl} = 1 | Y_{cl} = y),
   * more generally, a numerical entry for the given cell and locus.
   */
  def double get(Cell cell, Locus locus) 
  def Set<Cell> getCells() // assume stable iterator order
  def Set<Locus> getLoci()
  
  /**
   * In the context of constructing a likelihood, log(Z),
   * more generally, a numerical entry for the given cell and locus.
   */
  def double logNormalization() { 0.0 }
}