package corrupt.post

import corrupt.Cell
import corrupt.Locus
import java.util.Set
import blang.inits.Implementations

@Implementations(ReadOnlyCLMatrix, NoisyBinaryCLMatrix)
interface CellLocusMatrix {
  /**
   * Assuming an artificial uniform prior, return a  
   * "local" posterior probability that the given cell has the 
   * trait at the given locus. 
   */
  def double getTipAsDouble(Cell cell, Locus locus) 
  def Set<Cell> getCells() // assume stable iterator order
  def Set<Locus> getLoci()
}