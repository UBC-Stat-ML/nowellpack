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
import java.util.List
import corrupt.viz.MatrixViz.CellFiller

class PerfectPhyloViz extends Viz {
  val TreeViz<Set<TreeNode>> treeViz
  val MatrixViz matrixViz
  
  new (PerfectPhylo phylo, List<CellLocusMatrix> matrices) {
    if (matrices.size === 0 || matrices.size > 4)
      throw new RuntimeException
    val int groupSize = matrices.size + 1 // keep one to space out the groups
    val collapsedTree = phylo.collapsedTree 
    this.treeViz = new TreeViz(collapsedTree) => [ declareHeight(1.0) ]
    val CellLocusMatrix matrix = matrices.get(0) 
    val converted = MatrixOperations::dense(matrix.cells.size, matrix.loci.size * groupSize) 
    var int colIndex = 0
    for (locus : matrix.loci) {
      for (entry : treeViz.tipIndices.entrySet) {
        val cells = entry.key.filter[it instanceof Cell]
        if (cells.size != 1) 
          throw new RuntimeException
        val cell = cells.iterator.next as Cell 
        for (i : 0 ..< matrices.size)
          converted.set(entry.value, groupSize * colIndex + i, matrices.get(i).getTipAsDouble(cell, locus))
      }
      colIndex++
    }
    val CellFiller filler = [r, c, v, result | 
      val modulo = c % groupSize
      if (modulo === 0 || modulo === groupSize - 1)
        MatrixViz::greyScale.colour(0, 0, v, result)
      else 
        MatrixViz::colours(modulo - 1).colour(0, 0, v, result)
    ]
    this.matrixViz = new MatrixViz(converted, filler) => [ declareHeight(1.0) ]
  }
  
  override protected draw() {
    addChild(treeViz, 0, 0)
    addChild(matrixViz, treeViz.inferWidth, 0.0f)
  }
  
  override protected size() {
    return (treeViz.inferWidth + matrixViz.inferWidth) -> (treeViz.inferHeight)
  }
  
  public static def void main(String [] args) { 
    val rand = new Random(1)
    val phylo = new PerfectPhylo(CorruptStaticUtils::syntheticCells(100), CorruptStaticUtils::syntheticLoci(100))
    phylo.sampleUniform(new Random(1))
    val indicators = CLMatrixUtils::fromPhylo(phylo)
    val data = CLMatrixUtils::syntheticInclusionPrs(rand, phylo, 0.4)
    new PerfectPhyloViz(phylo, #[indicators, data]) => [
      declareHeight(800) 
      show
    ]
  }
}