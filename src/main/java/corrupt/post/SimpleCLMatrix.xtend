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
import bayonet.math.SpecialFunctions

@Data class SimpleCLMatrix implements CellLocusMatrix {
  
  val Indexer<Cell> cellsIdx
  val Indexer<Locus> lociIdx
  val Matrix matrix 
  
  new(Collection<Cell> cells, Collection<Locus> loci) {
    this.cellsIdx = new Indexer(GenomeMap::sanitize(cells))
    this.lociIdx = new Indexer(GenomeMap::orderLoci(loci))
    this.matrix = MatrixOperations::dense(cellsIdx.size, lociIdx.size)
  }
  
  override get(Cell cell, Locus locus)       { matrix.get(cellsIdx.o2i(cell), lociIdx.o2i(locus)) }
  def void set(Cell cell, Locus locus, double value) { 
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
        increment(cell, locus, another.get(cell, locus))
  }
  
  def void +=(PerfectPhylo phylo) {
    for (locus : loci) {
      val tips = phylo.getTips(locus)
      for (entry : tips.entrySet) 
        increment(entry.key, locus, if (entry.value) 1.0 else 0.0)
    }
  }
  
  def void increment(Cell cell, Locus locus, double increment) {
    set(cell, locus, get(cell, locus) + increment)
  }
  
  def void /=(double divisor) {
    for (cell : cells) 
      for (locus : loci) 
        set(cell, locus, get(cell, locus) / divisor) 
  }
  
  def logisticTransform() {
    for (var int r = 0; r < matrix.nRows; r++) 
      for (var int c = 0; c < matrix.nCols; c++) {
        val rho = matrix.get(r,c)
        matrix.set(r, c, SpecialFunctions::logistic(2*rho - 1))
      }
  }
  
}