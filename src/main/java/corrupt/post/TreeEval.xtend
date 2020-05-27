package corrupt.post

import corrupt.PerfectPhylo
import briefj.collections.Tree
import corrupt.Locus
import corrupt.Cell
import briefj.BriefIO
import java.io.File
import xlinear.Matrix

import static xlinear.MatrixOperations.*
import static extension xlinear.MatrixExtensions.*
import corrupt.post.CorruptPostProcessor.GoodnessOfFit
import java.util.Map
import java.util.LinkedHashMap
import java.util.Set

import static extension corrupt.post.CorruptPostProcessor.GoodnessOfFit.*
import blang.inits.experiments.Experiment
import blang.inits.Arg

class TreeEval extends Experiment {
  
  @Arg BinaryCLMatrix data
  
  @Arg File newickFile
    
  override run() {
    val treeStr = BriefIO::fileToString(newickFile)
    val sTree = PerfectPhylo::simpleTreeParsing(treeStr)
    val score = oraclePlacementScore(sTree, data)
    score.logTo(results.getTabularWriter("metrics")) 
    score.logToByLocus(results.getTabularWriter("metrics-by-locus"))
    results.getTabularWriter("tree-stats").write(
      "nNodes" -> sTree.postOrderTraversal.size,
      "depth" -> depth(sTree)
    )
  }
  
  def static <T> int depth(Tree<T> t) {
    if (t.children.empty) return 0
    return 1 + t.children.map[depth].max
  }
  
  def static Set<Cell> leaves(Tree<conifer.TreeNode> tree, Set<Cell> cells) {
    return tree.postOrderTraversal.map[new Cell(label.toString)].filter[cells.contains(it)].toSet
  }
  
  /**
   * For all possible clades (encoded by the node at the top), compute the score,
   * i.e. counts[truth][guess]
   * return map of clade -> counts
   * 
   * Not optimized for speed.
   */
  def static Map<conifer.TreeNode,Matrix> scores(Tree<conifer.TreeNode> simpleTree, Locus placed, BinaryCLMatrix data) {
    // check leaves match up 
    if (leaves(simpleTree, data.cells) != data.cells) {
      throw new RuntimeException
    }
    val result = new LinkedHashMap<conifer.TreeNode,Matrix>
    for (subTree : simpleTree.postOrderTraversal) {
      // create the clade
      val Set<Cell> clade = leaves(subTree, data.cells)
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
    scores(simpleTree, placed, data).values.maxBy[acc]
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
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}