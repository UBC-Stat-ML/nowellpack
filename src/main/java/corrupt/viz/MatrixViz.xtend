package corrupt.viz

import xlinear.Matrix
import org.eclipse.xtend.lib.annotations.Data
import xlinear.MatrixOperations
import processing.core.PApplet
import corrupt.viz.PublicSize
import corrupt.viz.Viz.PrivateSize

@Data class MatrixViz extends Viz  {
  val Matrix m
  val CellFiller filler
  
  new (Matrix m, CellFiller filler, PublicSize size) {
    super(size)
    this.m = m
    this.filler = filler
  }
  
  @FunctionalInterface
  static interface CellFiller {
    def void colour(int row, int col, double value, PApplet result)
  }
  
  // Assumes matrix entries between zero and one
  def static MatrixViz greyScale(Matrix m, PublicSize size) {
    return new MatrixViz(m, greyScale, size)
  }
  
  public static val CellFiller greyScale = [__, ___, v, result | result.fill(((1.0 - v) * 255.0) as int)]
  
  // Assumes matrix entries between zero and one
  static def CellFiller colours(int index, int paletteSize) {
    val float norm = paletteSize
    if (index < 0 || index >= paletteSize)
      throw new RuntimeException
    return [__, ___, v , result | 
      result.colorMode(PApplet::HSB, paletteSize) 
      result.fill(index as float, v as float * norm, norm)  
    ]
  }
  
  override draw() {
    noStroke
    pushStyle
    for (var int r = 0; r < m.nRows; r++) { 
      for (var int c = 0; c < m.nCols; c++) { 
        filler.colour(r, c, m.get(r, c), applet)
        rect(c, r, 1.0f, 1.0f)
      }
    }
    popStyle
  }
  
  override privateSize() { new PrivateSize(m.nCols as float,m.nRows as float) }
  
  public static def void main(String [] args) {
    val mtx = MatrixOperations::dense(1,10) 
    //val random = new Random(1)
    mtx.editInPlace[r, c, v| (r + c) / 10.0]
    for (i : 0 ..< 10)
      new MatrixViz(mtx, colours(i, 10), fixWidth(500)).output("" + i + ".pdf")   
  }
}