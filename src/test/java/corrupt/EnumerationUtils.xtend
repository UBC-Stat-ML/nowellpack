package corrupt

import java.util.List
import blang.runtime.SampledModel
import java.util.Set
import bayonet.distributions.ExhaustiveDebugRandom
import java.util.LinkedHashMap
import java.util.Map

import java.util.ArrayList
import static extension corrupt.CorruptExtensionUtils.*
import static corrupt.CorruptStaticUtils.*

class EnumerationUtils {
  
  def static List<SampledModel> enumerate(Set<Cell> cells, Set<Locus> loci) {
    val phylo = new PerfectPhylo(cells, loci)
    val current = new SampledModel(new Uniform.Builder().setPhylo(phylo).build)
    val result = new ArrayList
    val exhaustive = new ExhaustiveDebugRandom
    println("---")
    while (exhaustive.hasNext) {
      val copy = current.duplicate
      sampleNonUniform(exhaustive, (copy.model as Uniform).phylo)
      result.add(copy)
      println((current.model as Uniform).phylo)
    }
    println("---")
    return result 
  }
  
  def static private void sampleNonUniform(ExhaustiveDebugRandom rand, PerfectPhylo phylo) {
    for (locus : phylo.loci)
      phylo.tree.collapseEdge(locus)
    for (locus : phylo.loci) 
      splitSampler(phylo.tree, locus).sample(rand)
    phylo.updateAllSplits
  }
  
  def static private SplitSampler splitSampler(DirectedTree<TreeNode> phylogeny, Locus locusToAdd) {
    return new SplitSampler(phylogeny, locusToAdd, missingTips(phylogeny.cells)) 
  }
  
  def static private Map<Cell, SubtreeLikelihood> missingTips(List<Cell> cells) {
    val result = new LinkedHashMap
    for (cell : cells)
      result.put(cell, SubtreeLikelihood::missingTip)
    return result
  }
}