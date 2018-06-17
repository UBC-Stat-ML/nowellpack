package corrupt.viz

import corrupt.PerfectPhylo
import corrupt.post.CellLocusMatrix
import xlinear.MatrixOperations

import corrupt.CorruptStaticUtils
import bayonet.distributions.Random
import corrupt.post.CLMatrixUtils
import corrupt.TreeNode
import java.util.Set
import corrupt.Cell
import org.eclipse.xtext.resource.SaveOptions.Builder

class MatrixTreeViz extends Viz {
  val TreeViz<Set<TreeNode>> treeViz
  val MatrixViz matrixViz
  
  new (PerfectPhylo phylo, CellLocusMatrix matrix) {
    val collapsedTree = phylo.collapsedTree 
    this.treeViz = new TreeViz(collapsedTree) => [ declareHeight(1.0) ]
    val converted = MatrixOperations::dense(matrix.cells.size, matrix.loci.size) 
    var int colIndex = 0
    for (locus : matrix.loci) {
      for (entry : treeViz.tipIndices.entrySet) {
        val cells = entry.key.filter[it instanceof Cell]
        if (cells.size != 1) 
          throw new RuntimeException
        val cell = cells.iterator.next as Cell 
        converted.set(entry.value, colIndex, 1.0 - matrix.getTipAsDouble(cell, locus))
      }
      colIndex++
    }
    this.matrixViz = new MatrixViz(converted) => [ declareHeight(1.0) ]
  }
  
  override protected draw() {
    addChild(treeViz, 0, 0)
    addChild(matrixViz, treeViz.inferWidth, 0.0f)
  }
  
  override protected size() {
    return (treeViz.inferWidth + matrixViz.inferWidth) -> (treeViz.inferHeight)
  }
  
  public static def void main(String [] args) { 
    val phylo = new PerfectPhylo(CorruptStaticUtils::syntheticCells(1000), CorruptStaticUtils::syntheticLoci(50))
    phylo.sampleUniform(new Random(1))
    val matrix = CLMatrixUtils::fromPhylo(phylo)
    new MatrixTreeViz(phylo, matrix) => [
      declareHeight(800) 
      output("test.pdf")
    ]
  }
}