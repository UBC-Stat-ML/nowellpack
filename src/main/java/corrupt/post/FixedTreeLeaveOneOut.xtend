package corrupt.post

import corrupt.PerfectPhylo
import briefj.collections.Tree
import corrupt.Locus
import corrupt.Cell
import corrupt.CorruptStaticUtils
import java.util.ArrayList
import briefj.BriefIO
import java.io.File
import xlinear.Matrix

import static xlinear.MatrixOperations.*
import static extension xlinear.MatrixExtensions.*
import corrupt.post.CorruptPostProcessor.GoodnessOfFit
import conifer.UnrootedTreeUtils
import java.util.Map
import java.util.LinkedHashMap
import java.util.Set

import static extension corrupt.post.CorruptPostProcessor.GoodnessOfFit.*

/**
 * Inputs:
 *  - tree where leaves are cells
 *  - binary matrix
 * Output: quality score
 * 
 * For each locus, averaging
 *    For each internal node, maximizing
 *        Compute accuracy
 * 
 * For each FPR, FNR, maximizing on a grid
 * For each locus l, averaging
 *    For each cell c, averaging
 *        Create matrix
 *        For each internal edge, maximizing matrix-based score
 *        Check if prediction OK for the heldout entry
 * 
 */
class FixedTreeLeaveOneOut {
  
  /**
   * For all possible clades (encoded by the node at the top), compute the score,
   * i.e. counts[truth][guess]
   * 
   * Not optimized for speed, ie. 
   */
  def static Map<conifer.TreeNode,Matrix> scores(Tree<conifer.TreeNode> simpleTree, Locus placed, BinaryCLMatrix data) {
    val result = new LinkedHashMap<conifer.TreeNode,Matrix>
    for (subTree : simpleTree.postOrderTraversal) {
      // create the clade
      val Set<Cell> clade = subTree.postOrderTraversal.filter[label.toString.startsWith(Cell::PREFIX )].map[new Cell(label.toString)].toSet
      // score
      val Matrix score = dense(2,2)
      for (cell : data.cells) {
        val truth = data.get(cell, placed) as int
        val guess = if (clade.contains(cell)) 1 else 0
        score.increment(truth, guess, 1.0)
      }
      val key = subTree.label
      result.put(key, score)
    }
    return result
  }
  
  def static Matrix oraclePlacementScore(Tree<conifer.TreeNode> simpleTree, Locus placed, BinaryCLMatrix data) {
    scores(simpleTree, placed, data).values.maxBy[gof]
  }
  
  def static GoodnessOfFit oraclePlacementScore(Tree<conifer.TreeNode> simpleTree, BinaryCLMatrix data) {
    val result = new GoodnessOfFit
    for (locus : data.loci) {
      val score = oraclePlacementScore(simpleTree, locus, data)
      result.byLocus.put(locus, score)
      result.counts += score
    }
    return result
  }
  
  /// stuff below does not work since assumes eg can get leave
  
  //def static void score(Tree<conifer.TreeNode> _simpleTree, )
  
//  def static Matrix score(PerfectPhylo tree, Locus placed, BinaryCLMatrix data, Cell heldout) {
//    val result = dense(2,2)
//    val tips = tree.getTips(placed)
//    for (entry : tips.entrySet) 
//      if (entry.key != heldout) {
//      val truth = data.get(entry.key, placed) as int
//      val guess = if (entry.value) 1 else 0
//      result.increment(truth, guess, 1.0)
//    }
//    return result
//  }
  
//  def static oraclePlacementScore(Tree<conifer.TreeNode> simpleTree, BinaryCLMatrix data) {
//    val gof = new GoodnessOfFit
//    for (locus : data.loci) {
//      
//      // loop over locations
//      for (candidate : simpleTree.postOrderTraversal) {
//        val node = candidate.label
//        
//      }
//      
//    }
//  }
  
//  def static PerfectPhylo placeLocus(Tree<conifer.TreeNode> _simpleTree, Locus placed, conifer.TreeNode location) {
//    // ensure has root label on top
//    val Tree<conifer.TreeNode> simpleTree = 
//      if (_simpleTree.label.toString != CorruptStaticUtils::root.toString)
//        new Tree<conifer.TreeNode>(conifer.TreeNode.withLabel(CorruptStaticUtils::root.toString), new ArrayList => [add(_simpleTree)])
//      else 
//        _simpleTree
//    val result = new PerfectPhylo(simpleTree) [ conifer.TreeNode node |
//      val description = node.toString
//      if (description.startsWith(Cell::PREFIX )) return new Cell (description)
//      if (description == location.toString) return placed
//      if (description == CorruptStaticUtils::root.toString) return CorruptStaticUtils::root
//      return new Locus(node.toString)
//    ]
//    if (!result.loci.contains(placed)) throw new RuntimeException
//    return result
//  }
  
  def static void main(String [] args) {
    val file = new File("/Users/bouchard/experiments/corrupt-nextflow/results/all/2020-05-25-19-54-33-Zhb2qanN.exec/consensus.newick")
    val treeStr = BriefIO::fileToString(file)
    val sTree = PerfectPhylo::simpleTreeParsing(treeStr)
    
    val data = BinaryCLMatrix::create("/Users/bouchard/experiments/corrupt-nextflow/work/ae/fbba539aabb2824b56c5e191942eaa/filtered.csv")
    
    val gof = oraclePlacementScore(sTree, data)
    
    println(gof.counts.gof)
    
//    val removedLocus = treeFile.replaceAll("locus", "")
//    println(removedLocus)
//    val simpleTree = PerfectPhylo::simpleTreeParsing(removedLocus)
//    val placed = placeLocus(simpleTree, new Locus("TEST"), conifer.TreeNode.withLabel("_1_51500001_52000000"))
//    println(placed)
  }
}