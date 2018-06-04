package corrupt

import java.util.Map
import java.util.LinkedHashMap
import briefj.Indexer

import static extension corrupt.CorruptUtils.lociAndRoot
import bayonet.distributions.Random

class SplitSampler {
  val DirectedTree<TreeNode> phylogeny 
  val Map<TreeNode, SubtreeLikelihood> likelihoods = new LinkedHashMap
  val double rootExclLog
  val Indexer<TreeNode> rootLociIndexer = new Indexer
  
  new (DirectedTree<TreeNode> phylogeny, Locus locusToAdd, Map<Cell, SubtreeLikelihood> cellLikelihoods) {
    this.phylogeny = phylogeny
    rootLociIndexer.addAllToIndex(phylogeny.lociAndRoot) 
    if (rootLociIndexer.containsObject(locusToAdd))
      throw new RuntimeException
    likelihoods.putAll(cellLikelihoods) // (1)
    rootExclLog = computeLikelihood(phylogeny.root).exclusionLog
  }
  
  /**
   * Sample from the likelihood times a uniform prior.
   */
  def void sample(Random random) {
    
    
//    val double [] prs = newDoubleArrayOfSize(rootLociIndexer.size)
//    for (parentIndex : 0 ..< rootLociIndexer.size)
//      prs.set(parentIndex, )
  }
  
  def private SubtreeLikelihood computeLikelihood(TreeNode node) {
    if (node instanceof Cell) 
      return likelihoods.get(node) // see (1)
    var sumLogIncl = 0.0
    var sumLogExcl = 0.0
    for (child : phylogeny.children(node)) {
      var childLikelihood = computeLikelihood(child)
      sumLogIncl += childLikelihood.inclusionLog
      sumLogExcl += childLikelihood.exclusionLog
    }
    val result = new SubtreeLikelihood(sumLogIncl, sumLogExcl)
    likelihoods.put(node, result)
    return result
  }
}