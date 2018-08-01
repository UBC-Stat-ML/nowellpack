package corrupt

import bayonet.distributions.ExhaustiveDebugRandom
import bayonet.distributions.Random
import blang.mcmc.internals.SamplerBuilder
import blang.runtime.Observations
import blang.runtime.SampledModel
import blang.runtime.internals.objectgraph.GraphAnalysis
import corrupt.post.CLMatrixUtils
import java.util.ArrayList
import java.util.LinkedHashMap
import java.util.List
import java.util.Map

import static corrupt.CorruptStaticUtils.*

import static extension corrupt.CorruptExtensionUtils.*
import static extension corrupt.SplitSampler.sampleInPlace

class EnumerationUtils {
  
  def static syntheticModel(Random rand, int nCells, int nLoci) {
    val cells = syntheticCells(nCells)
    val loci = syntheticLoci(nLoci)
    val phylo = new PerfectPhylo(cells, loci)
    phylo.sampleUniform(rand)
    val tipInclPrs = CLMatrixUtils::syntheticInclusionPrs(rand, phylo, 0.5) 
    println(tipInclPrs)
    return new FixedMatrixModel.Builder()
      .setTipInclusionProbabilities(tipInclPrs)
      .build 
  }
  
  def static List<SampledModel> enumerateSyntheticModels(int nCells, int nLoci, double annealParam) {
    val rand = new Random(1)
    val model = syntheticModel(rand, nCells, nLoci)
    val analysis = new GraphAnalysis(model)
    val samplers = SamplerBuilder::build(analysis)
    var sModel = new SampledModel(analysis, samplers, true, true, null)
    sModel.exponent = annealParam 
    
    val result = new ArrayList
    val exhaustive = new ExhaustiveDebugRandom
    while (exhaustive.hasNext) {
      // freeze observation
      var copy = sModel.duplicate
      val observations = new Observations
      sampleNonUniform(exhaustive, (copy.model as FixedMatrixModel).phylo.getReconstruction)  
      (copy.model as FixedMatrixModel).phylo.cache.reset
      val graphAnalysis = new GraphAnalysis(copy.model, observations)
      copy = new SampledModel(graphAnalysis, SamplerBuilder.build(graphAnalysis), true, true, null)
      copy.exponent = annealParam 
      result.add(copy)
    }
    return result 
  }
  
  def static private void sampleNonUniform(ExhaustiveDebugRandom rand, PerfectPhylo phylo) {
    for (locus : phylo.loci)
      phylo.tree.collapseEdge(locus) 
    for (locus : phylo.loci) 
      phylo.tree.sampleInPlace(locus, missingTips(phylo.tree.cells), rand)
  }
  
  def static private Map<Cell, SubtreeLikelihood> missingTips(List<Cell> cells) {
    val result = new LinkedHashMap
    for (cell : cells)
      result.put(cell, SubtreeLikelihood::missingTip)
    return result
  }
}