package corrupt.post

import corrupt.Cell
import corrupt.Locus
import org.eclipse.xtend.lib.annotations.Data
import briefj.Indexer
import xlinear.Matrix
import xlinear.MatrixOperations
import java.util.Collection
import java.io.File
import briefj.BriefIO

@Data class SimpleCLMatrix implements CellLocusMatrix {
  
  val Indexer<Cell> cellsIdx
  val Indexer<Locus> lociIdx
  val Matrix matrix
  
  new(Collection<Cell> cells, Collection<Locus> loci) {
    this.cellsIdx = new Indexer(cells)
    this.lociIdx = new Indexer(loci)
    this.matrix = MatrixOperations::dense(cellsIdx.size, lociIdx.size)
  }
  
  override getTipAsDouble(Cell cell, Locus locus)       { matrix.get(cellsIdx.o2i(cell), lociIdx.o2i(locus)) }
  def void setTip(Cell cell, Locus locus, double value) { 
    matrix.set(cellsIdx.o2i(cell), lociIdx.o2i(locus), value)
  }
  
  override getCells() { cellsIdx.objects }
  override getLoci()  { lociIdx.objects}
  
  def void +=(CellLocusMatrix another) {
    CLMatrixUtils::checkCompatible(this, another)
    for (cell : cells) 
      for (locus : loci) 
        setTip(cell, locus, getTipAsDouble(cell, locus) + another.getTipAsDouble(cell, locus))
  }
  
  def void /=(double divisor) {
    for (cell : cells) 
      for (locus : loci) 
        setTip(cell, locus, getTipAsDouble(cell, locus) / divisor)
  }
}