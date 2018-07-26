package corrupt.post

import corrupt.Cell
import corrupt.Locus
import org.eclipse.xtend.lib.annotations.Data
import briefj.Indexer
import xlinear.Matrix
import xlinear.MatrixOperations
import java.util.Collection
import corrupt.PerfectPhylo
import corrupt.GenomeMap

@Data class SimpleCLMatrix implements CellLocusMatrix {
  
  val Indexer<Cell> cellsIdx
  val Indexer<Locus> lociIdx
  val Matrix matrix 
  
  new(Collection<Cell> cells, Collection<Locus> loci) {
    this.cellsIdx = new Indexer(GenomeMap::sanitize(cells))
    this.lociIdx = new Indexer(GenomeMap::orderLoci(loci))
    this.matrix = MatrixOperations::dense(cellsIdx.size, lociIdx.size)
  }
  
  override getTipAsDouble(Cell cell, Locus locus)       { matrix.get(cellsIdx.o2i(cell), lociIdx.o2i(locus)) }
  def void setTip(Cell cell, Locus locus, double value) { 
    matrix.set(cellsIdx.o2i(cell), lociIdx.o2i(locus), value)
  }
  
  def Matrix slice(Cell cell) {
    matrix.row(cellsIdx.o2i(cell))
  }
  
  def Matrix slice(Locus locus) {
    matrix.col(lociIdx.o2i(locus))
  }
  
  override getCells() { cellsIdx.objects }
  override getLoci()  { lociIdx.objects}
  
  def void +=(SimpleCLMatrix another) {
    CLMatrixUtils::checkCompatible(this, another)
    for (cell : cells) 
      for (locus : loci) 
        increment(cell, locus, another.getTipAsDouble(cell, locus))
  }
  
  def void +=(PerfectPhylo phylo) {
    for (locus : loci) {
      val tips = phylo.getTips(locus)
      for (entry : tips.entrySet) 
        increment(entry.key, locus, if (entry.value) 1.0 else 0.0)
    }
  }
  
  def void increment(Cell cell, Locus locus, double increment) {
    setTip(cell, locus, getTipAsDouble(cell, locus) + increment)
  }
  
  def void /=(double divisor) {
    for (cell : cells) 
      for (locus : loci) 
        setTip(cell, locus, getTipAsDouble(cell, locus) / divisor)
  }
}