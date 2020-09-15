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
  
  /**
   * Sample from the likelihood times a uniform prior, by doing block 
   * Gibbs sampling of one column of the cell-locus matrix.
   * Modifies the phylogeny in place.
   */
  def static SplitSampler sampleInPlace(
    DirectedTree<TreeNode> phylogeny, 
    Locus locusToAdd, 
    Map<Cell, SubtreeLikelihood> cellLikelihoods, Random rand
  ) {
    val result = new SplitSampler(phylogeny, cellLikelihoods)
    result.sample(rand, locusToAdd)
    if (result.rootLociIndexer.containsObject(locusToAdd))
      throw new RuntimeException
    return result
  }
  
  def static double maxLogConditional(
    DirectedTree<TreeNode> phylogeny, 
    Map<Cell, SubtreeLikelihood> cellLikelihoods
  ) {
    new SplitSampler(phylogeny, cellLikelihoods).maxLogConditional().value
  }
  
  /**
   * Return log probability of the argmax
   */
  def static double maximizeInPlace(
    DirectedTree<TreeNode> phylogeny, 
    Locus locusToAdd, 
    Map<Cell, SubtreeLikelihood> cellLikelihoods
  ) {
    val sampler = new SplitSampler(phylogeny, cellLikelihoods)
    if (sampler.rootLociIndexer.containsObject(locusToAdd))
      throw new RuntimeException
    return sampler.maximize(locusToAdd)
  }
  
  val DirectedTree<TreeNode> phylogeny 
  val Map<TreeNode, SubtreeLikelihood> likelihoods = new LinkedHashMap
  
  @Accessors(PUBLIC_GETTER)
  val Indexer<TreeNode> rootLociIndexer = new Indexer
  
  new (DirectedTree<TreeNode> phylogeny, Map<Cell, SubtreeLikelihood> cellLikelihoods) {
    this.phylogeny = phylogeny
    rootLociIndexer.addAllToIndex(phylogeny.lociAndRoot) 
    likelihoods.putAll(cellLikelihoods) // (1)
    computeLikelihood(phylogeny.root)
  }
  
  /**
   * Build the conditional probability for each possible attachment. 
   * Return the maximum and the index of the corresponding attachment.
   */
  private def Pair<Integer,Double> maxLogConditional() {
    val double [] logPrs = newDoubleArrayOfSize(rootLociIndexer.size)
    for (parentIndex : 0 ..< rootLociIndexer.size) {
      val node = rootLociIndexer.i2o(parentIndex)
      var attachPr = logAttachPr(node)
      val children = phylogeny.children(node)
      for (child : children) {
        val recursions = likelihoods.get(child)
        var inclPr = inclusionPr(recursions)
        inclPr = Math.max(inclPr, 1.0 - inclPr)
        attachPr += Math.log(inclPr)
      }
      logPrs.set(parentIndex, attachPr)
    }
    var double max = Double.NEGATIVE_INFINITY
    var int argmax = -1
    for (parentIndex : 0 ..< rootLociIndexer.size)
      if (logPrs.get(parentIndex) > max) {
        max = logPrs.get(parentIndex)
        argmax = parentIndex
      }
    return argmax -> max
  }
  
  /**
   * Return log probability of the argmax
   */
  private def double maximize(Locus locusToAdd) {
    val pair = maxLogConditional
    val pickedParent = rootLociIndexer.i2o(pair.key) 
    val children = phylogeny.children(pickedParent)
    val childrenToMove = new ArrayList
    for (child : children) {
      val recursions = likelihoods.get(child)
      if (inclusionPr(recursions) > 0.5) 
        childrenToMove.add(child)
    }
    phylogeny.addEdge(pickedParent, locusToAdd, childrenToMove)
    return pair.value
  }
  
  def double[] attachProbabilities() {
    val prs = newDoubleArrayOfSize(rootLociIndexer.size)
    for (parentIndex : 0 ..< rootLociIndexer.size) {
      val node = rootLociIndexer.i2o(parentIndex)
      prs.set(parentIndex, logAttachPr(node))
    }
    expNormalize(prs)
    return prs
  }
  
  private def void sample(Random random, Locus locusToAdd) {
    // step 1 of the sampling process: pick a parent    
    val double [] prs = attachProbabilities()
    val sampledParent = rootLociIndexer.i2o(random.nextCategorical(prs)) 
    // step 2 of the sampling process: which children to move?
    val children = phylogeny.children(sampledParent)
    val childrenToMove = new ArrayList
    for (child : children) {
      val recursions = likelihoods.get(child)
      val include = random.nextBernoulli(inclusionPr(recursions))
      if (include)
        childrenToMove.add(child)
    }
    // step 3: apply the change to the tree
    phylogeny.addEdge(sampledParent, locusToAdd, childrenToMove)
  }
  
  /**
   * This computes, in the paper's notation, "\bar \rho_v" up to a 
   * normalization constant.
   * 
   * Give the result in log scale.
   * 
   * The input of the function, "node" is "v" in the paper's 
   * mathematical notation.
   */  
  def private double logAttachPr(TreeNode node) {
    val recursions = likelihoods.get(node)
    return recursions.logProductPQ - recursions.exclusionLog
  }
  
  /**
   * Last equation of section 4.1 in latex document.
   */
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