package corrupt.viz

import blang.inits.ConstructorArg
import blang.inits.DesignatedConstructor
import corrupt.Cell
import corrupt.PerfectPhylo
import corrupt.TreeNode
import corrupt.post.CellLocusMatrix
import corrupt.viz.MatrixViz.CellFiller
import java.util.List
import java.util.Set
import xlinear.MatrixOperations
import blang.inits.experiments.Experiment
import corrupt.post.CLMatrixUtils

class PerfectPhyloViz extends Viz {
  val TreeViz<Set<TreeNode>> treeViz
  val MatrixViz matrixViz
  
  @DesignatedConstructor
  new (@ConstructorArg("phylo") PerfectPhylo phylo, @ConstructorArg("matrices") List<CellLocusMatrix> matrices, @ConstructorArg("size") PublicSize size) {
    super(size)
    
    if (matrices.size > 4)
      throw new RuntimeException
    
    // add the indicators automatically
    val indicators = CLMatrixUtils::fromPhylo(phylo)
    matrices.add(0, indicators)
    
    // setup tree  
    val collapsedTree = phylo.collapsedTree 
    this.treeViz = new TreeViz(collapsedTree, fixHeight(1))
     
    // overlay matrices 
    val int groupSize = matrices.size + 1 // keep one to space out the groups
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
      if (modulo < 2 || modulo === groupSize - 1)
        MatrixViz::greyScale.colour(0, 0, v, result)
      else 
        MatrixViz::colours(modulo - 2).colour(0, 0, v, result)
    ]
    this.matrixViz = new MatrixViz(converted, filler, fixHeight(1)) 
  }
  
  override protected draw() {
    addChild(treeViz, 0, 0)
    addChild(matrixViz, treeViz.publicWidth, 0.0f)
  }
  
  override protected privateSize() {
    return new PrivateSize(treeViz.publicWidth + matrixViz.publicWidth, treeViz.publicHeight)
  }
  
  public static def void main(String [] args) { 
    val code = Experiment.start(args) 
    if (code !== 0)
      System.exit(code)
    else
      {} // Do not exit; drawing thread probably still working
  }
}