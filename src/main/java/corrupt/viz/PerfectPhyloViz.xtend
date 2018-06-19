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
import java.util.Optional
import java.util.ArrayList

class PerfectPhyloViz extends Viz {
  val TreeViz<Set<TreeNode>> treeViz
  val MatrixViz matrixViz
  
  new (
    PerfectPhylo phylo, 
    List<CellLocusMatrix> matrices, 
    PublicSize size) {
    this(phylo, matrices, size, Optional.empty, Optional.empty)
  }
  
  new (
    PerfectPhylo phylo, 
    List<CellLocusMatrix> matrices, 
    PublicSize size,
    PerfectPhylo refPhylo
    ) {
    this(phylo, matrices, size, Optional.of(refPhylo), Optional.empty)
  }
  
  @DesignatedConstructor
  new (
    @ConstructorArg("phylo") PerfectPhylo phylo, 
    @ConstructorArg("matrices") List<CellLocusMatrix> matrices, 
    @ConstructorArg("size") PublicSize size,
    @ConstructorArg("ref") Optional<PerfectPhylo> refPhylo,
    @ConstructorArg(value = "colourCodes", description = coloursDescriptions) Optional<List<Integer>> codes
  ) {
    super(size)
    
    // add the indicators for the displayed tree + one optional reference tree
    val indicators = CLMatrixUtils::fromPhylo(phylo)
    matrices.add(0, indicators)
    
    if (refPhylo.present) {
      val refIndics = CLMatrixUtils::fromPhylo(refPhylo.get)
      matrices.add(1, refIndics)
    }
    
    // prepare colour schemes for matrices
    val schemes = schemes(matrices, codes)
    
    // setup tree  
    val collapsedTree = phylo.collapsedTree 
    this.treeViz = new TreeViz(collapsedTree, fixHeight(1))
     
    // overlay matrices 
    val int groupSize = schemes.size
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
    val CellFiller filler = [r, c, v, result | schemes.get(c % schemes.size).colour(0, 0, v, result) ]
    this.matrixViz = new MatrixViz(converted, filler, fixHeight(1)) 
  }
  
  private static def List<CellFiller> schemes(List<CellLocusMatrix> matrices, Optional<List<Integer>> _specifiedOrder) {
    val result = new ArrayList
    if (_specifiedOrder.present) {
      val specifiedOrder = _specifiedOrder.get
      if (specifiedOrder.size !== matrices.size) 
        throw new RuntimeException("The provided colour schemes should be of length " + matrices.size)
      for (i : specifiedOrder)  
        result.add(colourCodes.get(i))
    } else {
      for (i : 0..<matrices.size)
        result.add(colourCodes.get(i % colourCodes.size))
    }
    result.add(MatrixViz::greyScale) // hack for the white divisor between groups
    return result
  }
  
  static val coloursDescriptions = "list of integers where 0 is greyScale and 1--10 are colour palettes"
  public static val List<CellFiller> colourCodes = new ArrayList => [
    add(MatrixViz::greyScale)
    val nColours = 10
    for (i : 0 ..< nColours)
      add(MatrixViz::colours(i, nColours))
  ]
  
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