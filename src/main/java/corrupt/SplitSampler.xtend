package corrupt

import java.util.Map
import java.util.LinkedHashMap
import briefj.Indexer

import bayonet.distributions.Random
import static java.lang.Math.exp

import static bayonet.distributions.Multinomial.expNormalize 

import static bayonet.math.NumericalUtils.logAdd
import java.util.ArrayList
import static extension corrupt.CorruptExtensionUtils.*
import org.eclipse.xtend.lib.annotations.Accessors

class SplitSampler {
  
  @Accessors(PUBLIC_GETTER) var double logNormalization
  
  /**
   * Sample from the likelihood times a uniform prior. 
   * Modifies the phylogeny in place.
   */
  def static SplitSampler sampleInPlace(
    DirectedTree<TreeNode> phylogeny, 
    Locus locusToAdd, 
    Map<Cell, SubtreeLikelihood> cellLikelihoods, Random rand
  ) {
    val result = new SplitSampler(phylogeny, locusToAdd, cellLikelihoods)
    result.sample(rand)
    return result
  }
  
  val DirectedTree<TreeNode> phylogeny 
  val Map<TreeNode, SubtreeLikelihood> likelihoods = new LinkedHashMap
  val double rootExclLog
  val Indexer<TreeNode> rootLociIndexer = new Indexer
  val Locus locusToAdd
  
  private new (DirectedTree<TreeNode> phylogeny, Locus locusToAdd, Map<Cell, SubtreeLikelihood> cellLikelihoods) {
    this.locusToAdd = locusToAdd
    this.phylogeny = phylogeny
    rootLociIndexer.addAllToIndex(phylogeny.lociAndRoot) 
    if (rootLociIndexer.containsObject(locusToAdd))
      throw new RuntimeException
    likelihoods.putAll(cellLikelihoods) // (1)
    rootExclLog = computeLikelihood(phylogeny.root).exclusionLog
  }
  
  private def void sample(Random random) {
    // parent    
    val double [] prs = newDoubleArrayOfSize(rootLociIndexer.size)
    for (parentIndex : 0 ..< rootLociIndexer.size) {
      val node = rootLociIndexer.i2o(parentIndex)
      val recursions = likelihoods.get(node)
      val likelihood = rootExclLog + recursions.logProductPQ - recursions.exclusionLog
      prs.set(parentIndex, likelihood)
    }
    logNormalization = expNormalize(prs)
    val sampledParent = rootLociIndexer.i2o(random.nextCategorical(prs)) 
    // children to move
    val children = phylogeny.children(sampledParent)
    val childrenToMove = new ArrayList
    for (child : children) {
      val recursions = likelihoods.get(child)
      val include = random.nextBernoulli(inclusionPr(recursions))
      if (include)
        childrenToMove.add(child)
    }
    phylogeny.addEdge(sampledParent, locusToAdd, childrenToMove)
  }
  
  def private double inclusionPr(SubtreeLikelihood recursion) {
    return exp(recursion.inclusionLog - logAdd(recursion.inclusionLog, recursion.exclusionLog))
  }
  
  def private SubtreeLikelihood computeLikelihood(TreeNode node) {
    if (node instanceof Cell) 
      return likelihoods.get(node) // see (1)
    var sumLogIncl = 0.0
    var sumLogExcl = 0.0
    var sumLogChildrenPQ = 0.0
    for (child : phylogeny.children(node)) {
      var childLikelihood = computeLikelihood(child)
      sumLogIncl += childLikelihood.inclusionLog
      sumLogExcl += childLikelihood.exclusionLog
      sumLogChildrenPQ += logAdd(childLikelihood.inclusionLog, childLikelihood.exclusionLog)
    }
    val result = new SubtreeLikelihood(sumLogIncl, sumLogExcl, sumLogChildrenPQ)
    likelihoods.put(node, result)
    return result
  }
}