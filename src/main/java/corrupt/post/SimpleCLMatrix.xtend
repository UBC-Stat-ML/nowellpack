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
import corrupt.PerfectPhylo
import java.util.LinkedHashSet
import java.util.ArrayList
import java.util.Collections
import java.util.Comparator
import java.util.List

@Data class SimpleCLMatrix implements CellLocusMatrix {
  
  val Indexer<Cell> cellsIdx
  val Indexer<Locus> lociIdx
  val Matrix matrix
  
  new(Collection<Cell> cells, Collection<Locus> loci) {
    this.cellsIdx = new Indexer(sanitize(cells))
    this.lociIdx = new Indexer(sanitize(loci))
    this.matrix = MatrixOperations::dense(cellsIdx.size, lociIdx.size)
  }
  
  private static def  <T> List<T> sanitize(Collection<T> items) {
    val asSet = new LinkedHashSet(items)
    val asList = new ArrayList(asSet)
    Collections::sort(asList, Comparator.comparing[it.toString])
    return asList  
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