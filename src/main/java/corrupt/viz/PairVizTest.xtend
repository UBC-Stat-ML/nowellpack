package corrupt.viz

import xlinear.Matrix
import org.eclipse.xtend.lib.annotations.Data
import xlinear.MatrixOperations
import java.util.Random

@Data class PairVizTest extends Viz {
  
  val Matrix m1
  val Matrix m2
  
  override protected draw() {
    addChild(new MatrixViz(m1) => [declareWidth(1.0f)], 0.0f, 0.0f)
    addChild(new MatrixViz(m2) => [declareWidth(1.0f)], 1.0f, 1.0f)
  }
  
  override protected size() {
    Pair.of(2.0f, 2.0f)
  }
  
  public static def void main(String [] args) {
    val mtx = MatrixOperations::dense(20,30)
    val random = new Random(1)
    mtx.editInPlace[r, c, v| random.nextDouble]
    
    new PairVizTest(mtx, mtx * 2.0).show
  }
}