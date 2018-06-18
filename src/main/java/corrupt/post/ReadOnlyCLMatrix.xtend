package corrupt.post

import com.rits.cloning.Immutable
import org.eclipse.xtend.lib.annotations.Data
import corrupt.Cell
import corrupt.Locus

@Immutable
@Data class ReadOnlyCLMatrix implements CellLocusMatrix {
  val CellLocusMatrix enclosed 
  
  private new(CellLocusMatrix enclosed) {
    this.enclosed = enclosed
  }
  
  def static ReadOnlyCLMatrix readOnly(CellLocusMatrix m) {
    if (m instanceof ReadOnlyCLMatrix) return m
    else return new ReadOnlyCLMatrix(m)
  }
  
  override getTipAsDouble(Cell cell, Locus locus) { enclosed.getTipAsDouble(cell, locus) }
  override getCells() { enclosed.cells }
  override getLoci() { enclosed.loci}
}