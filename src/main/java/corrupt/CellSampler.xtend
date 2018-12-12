package corrupt

import bayonet.distributions.Random
import java.util.Map
import briefj.Indexer

import static bayonet.distributions.Multinomial.expNormalize 

import static extension corrupt.CorruptExtensionUtils.*
import static corrupt.CorruptStaticUtils.*
import java.util.Collections

class CellSampler {
  
  /**
   * Sample from the likelihood times a uniform prior, by doing block 
   * Gibbs sampling of one row of the cell-locus matrix.
   * Modifies the phylogeny in place.
   */
  def static CellSampler sampleInPlace(
    DirectedTree<TreeNode> phylogeny, 
    Cell cellToAdd,
    Map<Locus, SubtreeLikelihood> cellLikelihoods, 
    Random rand
  ) {
    val result = new CellSampler(phylogeny, cellLikelihoods)
    result.sample(rand, cellToAdd)
    return result
  }
  
  def static CellSampler maximizeInPlace(
    DirectedTree<TreeNode> phylogeny, 
    Cell cellToAdd,
    Map<Locus, SubtreeLikelihood> cellLikelihoods 
  ) {
    val result = new CellSampler(phylogeny, cellLikelihoods)
    result.maximize(cellToAdd) 
    return result
  }
  
  val DirectedTree<TreeNode> phylogeny 
  val Map<Locus, SubtreeLikelihood> likelihoods
  val Indexer<TreeNode> rootLociIndexer = new Indexer
  
  private new (DirectedTree<TreeNode> phylogeny, Map<Locus, SubtreeLikelihood> likelihoods) {
    this.phylogeny = phylogeny
    this.likelihoods = likelihoods
    rootLociIndexer.addAllToIndex(phylogeny.lociAndRoot) 
  }
  
  def private void sample(Random random, Cell cellToAdd) {
    val double [] prs = newDoubleArrayOfSize(rootLociIndexer.size)
    computeLogPrs(phylogeny.root, 0.0, prs)
    expNormalize(prs)
    val sampledParent = rootLociIndexer.i2o(random.nextCategorical(prs)) 
    phylogeny.addEdge(sampledParent, cellToAdd, Collections::emptyList)
  }
  
  def private void maximize(Cell cellToAdd) {
    val double [] prs = newDoubleArrayOfSize(rootLociIndexer.size)
    computeLogPrs(phylogeny.root, 0.0, prs)
    var max = Double::NEGATIVE_INFINITY
    var int argmax = -1
    for (i : 0 ..< prs.length) {
      val cur = prs.get(i)
      if (cur > max) {
        max = cur
        argmax = i
      }
    }
    val parent = rootLociIndexer.i2o(argmax) 
    phylogeny.addEdge(parent, cellToAdd, Collections::emptyList)
  }
  
  def private void computeLogPrs(TreeNode node, double parentLogPr, double [] prs) {
    if (node instanceof Cell) return
    val index = rootLociIndexer.o2i(node)
    val currentLogPr = parentLogPr + 
      if (node == root) 
        0.0
      else {
        val ll = likelihoods.get(node)
        ll.inclusionLog - ll.exclusionLog
      }
    prs.set(index, currentLogPr)
    for (TreeNode child : phylogeny.children(node))
      computeLogPrs(child, currentLogPr, prs)
  }
  
}

