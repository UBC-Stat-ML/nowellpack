package corrupt.post

import corrupt.Cell
import corrupt.Locus
import briefj.Indexer
import java.util.BitSet
import com.rits.cloning.Immutable
import blang.inits.DesignatedConstructor
import blang.inits.Input
import java.io.File

@Immutable
class BinaryCLMatrix implements CellLocusMatrix {
  val Indexer<Cell> cellsIdx
  val Indexer<Locus> lociIdx
  val BitSet bits
  
  def static BinaryCLMatrix binary(CellLocusMatrix m) {
    if (m instanceof BinaryCLMatrix) return m
    else return new BinaryCLMatrix(m)
  }
  
  private new (CellLocusMatrix model) {
    this.cellsIdx = new Indexer(model.cells)
    this.lociIdx = new Indexer(model.loci)
    this.bits = new BitSet(cellsIdx.size * lociIdx.size)
    for (cell : model.cells)
      for (locus : model.loci)
        switch (model.get(cell, locus)) {
          case 1.0 : bits.set(index(cell, locus))
          case 0.0 : {} // nothing to do
          default : throw new RuntimeException("Unexpected non binary entry: " + model.get(cell, locus))
        }
  }
  
  override get(Cell cell, Locus locus) {
    return if (bits.get(index(cell, locus))) 1.0 else 0.0
  }
  
  private def index(Cell cell, Locus locus) {
    cellsIdx.o2i(cell) + cellsIdx.size * lociIdx.o2i(locus)
  }
  
  override getCells() { cellsIdx.objects }
  override getLoci()  { lociIdx.objects}
  
  @DesignatedConstructor
  public static def BinaryCLMatrix create(
      @Input String path
  ) { 
    return new BinaryCLMatrix(CLMatrixUtils::fromCSV(new File(path)))    
  } 
}
