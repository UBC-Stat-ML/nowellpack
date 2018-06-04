package corrupt

import blang.mcmc.Sampler
import bayonet.distributions.Random
import blang.mcmc.SampledVariable
import blang.mcmc.internals.SamplerBuilderContext
import blang.mcmc.ConnectedFactor
import blang.core.LogScaleFactor
import java.util.List
import static blang.runtime.internals.objectgraph.StaticUtils.node
import java.util.Map
import java.util.LinkedHashMap

import static extension corrupt.CorruptUtils.cells

class PerfectPhyloGibbsSampler implements Sampler {
  @SampledVariable public Split split
  @ConnectedFactor public List<LogScaleFactor> _numericFactors
  Map<Cell, List<LogScaleFactor>> numericFactors
  
  override execute(Random rand) {
    split.tree.collapseEdge(split.locus)
    new SplitSampler(split.tree, split.locus, cellInclusionLogProbabilities).sample(rand)
    split.updateTips 
  }
  
  def Map<Cell,SubtreeLikelihood> cellInclusionLogProbabilities() {
    val result = new LinkedHashMap
    for (cell : split.tree.cells) {
      val logP = unnormalizedCellInclusionExclusionLogPr(cell, true)
      val logQ = unnormalizedCellInclusionExclusionLogPr(cell, false)
      result.put(cell, new SubtreeLikelihood(logP, logQ))
    }
    return result
  }
  
  def private double unnormalizedCellInclusionExclusionLogPr(Cell cell, boolean include) {
    val indic = split.tipIndicators.get(cell)
    val backup = indic.included
    indic.included = include
    val factors = numericFactors.get(cell)
    var result = 0.0
    for (factor : factors)
      result += factor.logDensity
    indic.included = backup
    return result
  }
  
  override setup(SamplerBuilderContext context) {
    // TODO: check for loops
    numericFactors = new LinkedHashMap
    for (indicator : split.tipIndicators.values) {
      val factors = context.connectedFactors(node(indicator))
      numericFactors.put(indicator.cell, factors.toArray as List<LogScaleFactor>)
    }
    return true
  }
}
