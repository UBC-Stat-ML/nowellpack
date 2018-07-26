package corrupt.post

import corrupt.Cell
import corrupt.Locus
import blang.inits.DesignatedConstructor
import blang.inits.Input
import java.io.File
import java.util.Set

interface CellLocusMatrix {
  def double getTipAsDouble(Cell cell, Locus locus) 
  def Set<Cell> getCells() // assume stable iterator order
  def Set<Locus> getLoci()
  
  @DesignatedConstructor
  public static def ReadOnlyCLMatrix create(
      @Input String path
  ) { 
    println("Loading tip inclusion probabilities")
    return ReadOnlyCLMatrix::readOnly(CLMatrixUtils::fromCSV(new File(path)))    
  } 
}