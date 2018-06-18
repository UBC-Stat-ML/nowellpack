package corrupt.viz

import xlinear.Matrix
import org.eclipse.xtend.lib.annotations.Data
import bayonet.distributions.Random
import xlinear.MatrixOperations
import processing.core.PApplet

@Data class MatrixViz extends Viz  {
  val Matrix m
  val CellFiller filler
  
  @FunctionalInterface
  static interface CellFiller {
    def void colour(int row, int col, double value, PApplet result)
  }
  
  // Assumes matrix entries between zero and one
  def static MatrixViz greyScale(Matrix m) {
    return new MatrixViz(m, greyScale)
  }
  
  public static val CellFiller greyScale = [__, ___, v, result | result.fill(((1.0 - v) * 255.0) as int)]
  
  // Assumes matrix entries between zero and one
  static def CellFiller colours(int index) {
    return [__, ___, v , result |  
      if (index < 0 || index > 2)
        throw new RuntimeException
      val varying = ((1.0 - v) * 255.0) as float
      val r = if (index === 0) 0f else varying
      val g = if (index === 1) 0f else varying
      val b = if (index === 2) 0f else varying
      result.fill(r, g, b)  
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
  
  override size() { Pair.of(m.nCols as float,m.nRows as float) }
  
  public static def void main(String [] args) {
    val mtx = MatrixOperations::dense(20,30) 
    val random = new Random(1)
    mtx.editInPlace[r, c, v| random.nextDouble]
    new MatrixViz(mtx, greyScale) => [
      declareWidth(1000)
      show
    ]
  }
}