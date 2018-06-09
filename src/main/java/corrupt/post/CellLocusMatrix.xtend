package corrupt.post

import corrupt.Cell
import corrupt.Locus
import java.util.Set

interface CellLocusMatrix {
  def double getTipAsDouble(Cell cell, Locus locus) 
  def Set<Cell> getCells()
  def Set<Locus> getLoci()
}