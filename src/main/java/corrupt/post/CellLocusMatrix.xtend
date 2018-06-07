package corrupt.post

import java.util.Collection
import corrupt.Cell
import corrupt.Locus

interface CellLocusMatrix {
  def double getTipAsDouble(Cell cell, Locus locus) 
  def Collection<Cell> getCells()
  def Collection<Locus> getLoci()
}