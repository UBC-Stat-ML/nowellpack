package corrupt

import java.util.List
import blang.runtime.SampledModel
import bayonet.distributions.ExhaustiveDebugRandom
import java.util.LinkedHashMap
import java.util.Map
import static extension corrupt.SplitSampler.sampleInPlace

import java.util.ArrayList
import static extension corrupt.CorruptExtensionUtils.*
import static corrupt.CorruptStaticUtils.*
import blang.io.GlobalDataSource
import blang.types.Plated
import blang.types.StaticUtils
import bayonet.distributions.Random
import blang.runtime.internals.objectgraph.GraphAnalysis
import blang.runtime.Observations
import blang.types.internals.SimplePlate
import blang.mcmc.internals.SamplerBuilder
import blang.types.Plate
import blang.core.RealVar
import corrupt.post.CLMatrixUtils

class EnumerationUtils {
  
  def static syntheticModel(Random rand, int nCells, int nLoci) {
    val cells = syntheticCells(nCells)
    val loci = syntheticLoci(nLoci)
    val phylo = new PerfectPhylo(cells, loci)
    phylo.sampleUniform(rand)
    val tipInclPrs = CLMatrixUtils::syntheticInclusionPrs(rand, phylo, 0.5) 
    println(tipInclPrs)
    return new CorruptModel.Builder()
      .setTipInclusionProbabilities(tipInclPrs)
      .build 
  }
  
  def static List<SampledModel> enumerateSyntheticModels(int nCells, int nLoci, double annealParam) {
    val rand = new Random(1)
    val model = syntheticModel(rand, nCells, nLoci)
    var sModel = new SampledModel(model)
    sModel.exponent = annealParam 
    
    val result = new ArrayList
    val exhaustive = new ExhaustiveDebugRandom
    while (exhaustive.hasNext) {
      // freeze observation
      var copy = sModel.duplicate
      val observations = new Observations
      sampleNonUniform(exhaustive, (copy.model as CorruptModel).phylo.phylo)  
      val graphAnalysis = new GraphAnalysis(copy.model, observations)
      copy = new SampledModel(graphAnalysis, SamplerBuilder.build(graphAnalysis), null)
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