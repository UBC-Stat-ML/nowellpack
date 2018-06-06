package corrupt.viz

import xlinear.Matrix
import org.eclipse.xtend.lib.annotations.Data
import bayonet.distributions.Random
import xlinear.MatrixOperations

@Data class MatrixViz extends Viz  {
  val Matrix m
  
  override draw() {
    noStroke
    for (var float r = 0; r < m.nRows; r++) { 
      for (var float c = 0; c < m.nCols; c++) { 
        fill((m.get(r as int, c as int) * 255.0) as int)
        rect(c, r, 1.0f, 1.0f)
      }
    }
  }
  
  override size() { Pair.of(m.nCols as float,m.nRows as float) }
  
  public static def void main(String [] args) {
    val mtx = MatrixOperations::dense(20,30)
    val random = new Random(1)
    mtx.editInPlace[r, c, v| random.nextDouble]
    new MatrixViz(mtx).show
  }
}