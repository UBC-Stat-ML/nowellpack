package corrupt.viz

import blang.inits.ConstructorArg
import blang.inits.DesignatedConstructor
import corrupt.Cell
import corrupt.PerfectPhylo
import corrupt.TreeNode
import corrupt.post.CellLocusMatrix
import viz.components.MatrixViz.CellFiller
import java.util.List
import java.util.Set
import xlinear.MatrixOperations
import blang.inits.experiments.Experiment
import corrupt.post.CLMatrixUtils
import java.util.Optional
import java.util.ArrayList
import viz.core.PublicSize
import viz.core.Viz
import viz.components.TreeViz
import viz.components.MatrixViz
import corrupt.GenomeMap
import java.util.Collection
import corrupt.Locus
import processing.core.PApplet
import corrupt.post.ReadOnlyCLMatrix
import java.util.LinkedHashSet
import java.io.File

class PerfectPhyloViz extends Viz {
  val TreeViz<Set<TreeNode>> treeViz
  val MatrixViz matrixViz
  val int nMatrices
  val Collection<Locus> loci
  
  new (
    PerfectPhylo phylo, 
    List<ReadOnlyCLMatrix> matrices, 
    PublicSize size) {
    this(phylo, matrices, size, Optional.empty, Optional.empty, Optional.empty)
  }
  
  new (
    PerfectPhylo phylo, 
    List<ReadOnlyCLMatrix> matrices, 
    PublicSize size,
    PerfectPhylo refPhylo
    ) {
    this(phylo, matrices, size, Optional.of(refPhylo), Optional.empty, Optional.empty)
  }
  
  def static void visualizePerChromosome(File directory, PerfectPhylo phylo, List<ReadOnlyCLMatrix> matrices, PublicSize size) {
    val allLoci = allLoci(matrices)
    val map = new GenomeMap(allLoci)
    directory.mkdirs
    for (chr : map.orderedChromosomes) {
      val viz = new PerfectPhyloViz(phylo, matrices, size, Optional.empty, Optional.empty, Optional.of(map.orderedLoci(chr).toSet))
      val outFile = new File(directory, chr + ".pdf") 
      viz.output(outFile)
    }
  } 
  
  //public static 
  
  def static List<Locus> allLoci(List<? extends CellLocusMatrix> matrices) {
    val loci = new LinkedHashSet<Locus>
    for (matrix : matrices)
      loci.addAll(matrix.loci)
    return GenomeMap::orderedLoci(loci)
  }
  
  @DesignatedConstructor
  new (
    @ConstructorArg("phylo") PerfectPhylo phylo, 
    @ConstructorArg("matrices") List<ReadOnlyCLMatrix> _matrices, 
    @ConstructorArg("size") PublicSize size,
    @ConstructorArg("ref") Optional<PerfectPhylo> refPhylo,
    @ConstructorArg(value = "colourCodes", description = coloursDescriptions) Optional<List<Integer>> codes,
    @ConstructorArg("restriction") Optional<Set<Locus>> restriction
  ) {
    super(size)
    
    // add the indicators for the displayed tree + one optional reference tree
    val indicators = CLMatrixUtils::fromPhylo(phylo)
    val ArrayList<CellLocusMatrix> matrices = new ArrayList(_matrices)
    matrices.add(0, indicators)
    
    if (refPhylo.present) {
      val refIndics = CLMatrixUtils::fromPhylo(refPhylo.get)
      matrices.add(1, refIndics)
    }
    
    // prepare colour schemes for matrices
    val schemes = schemes(matrices, codes)
    
    // setup tree  
    val collapsedTree = phylo.collapsedTree  
    this.treeViz = new TreeViz(collapsedTree.root, [collapsedTree.children(it)], fixHeight(1))  
     
    // overlay matrices 
    nMatrices = matrices.size
    val int groupSize = nMatrices + 1
    val allCells = matrices.get(0).cells
     
    loci = allLoci(matrices)
    if (restriction.present)
      loci.retainAll(restriction.get)
      
    val converted = MatrixOperations::dense(allCells.size, loci.size * groupSize) 
    var int colIndex = 0
    for (locus : loci) {
      for (entry : treeViz.tipIndices.entrySet) {
        val cells = entry.key.filter[it instanceof Cell]
        if (cells.size != 1) 
          throw new RuntimeException
        val cell = cells.iterator.next as Cell 
        for (i : 0 ..< matrices.size) {
          val value = if (matrices.get(i).loci.contains(locus)) matrices.get(i).get(cell, locus) else 0.0
          converted.set(entry.value, groupSize * colIndex + i, value)
        }
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
    
    // loci and chromosome markers
    strokeWeight(0.001f)
    for (var int i = 0; i < loci.size; i++)
      markGroup(i, i+1, 0.0f)
     
    pushStyle
    colorMode(PApplet::HSB, 25f) 
    if (GenomeMap::genomeMapFormatted(loci)) {
      val GenomeMap map = new GenomeMap(loci)
      var i = 0
      var c = 0
      for (chr : map.orderedChromosomes) {
        val nLoci = map.orderedLoci(chr).size
        stroke(c, 25.0f, 25.0f)
        markGroup(i, i+nLoci, 0.5f)
        i += nLoci
        c++
      }
    }
    popStyle
  }
  
  /**
   * Height between 0 and 1 relative to bottom
   */
  private def markGroup(int firstIncl, int lastExcl, float height) {
    val colWidth = matrixViz.publicWidth / matrixViz.m.nCols
    val y = treeViz.publicHeight + height * chrMarkerHeight
    val x1 = treeViz.publicWidth + firstIncl * colWidth * (nMatrices + 1)
    val x2 = treeViz.publicWidth + lastExcl * colWidth * (nMatrices + 1) - colWidth
    line(x1, y, x2, y)
  }
  static val chrMarkerHeight = 0.1f /* add 10% at bottom for chromosome markers */
  
  override protected privateSize() {
    return new PrivateSize(treeViz.publicWidth + matrixViz.publicWidth, treeViz.publicHeight + chrMarkerHeight) 
  }
  
  static def void main(String [] args) { 
    Experiment.startAutoExit(args)
  }
}