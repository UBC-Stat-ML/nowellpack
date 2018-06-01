package corrupt

import blang.mcmc.Sampler
import bayonet.distributions.Random
import blang.mcmc.SampledVariable
import blang.mcmc.internals.SamplerBuilderContext
import blang.mcmc.ConnectedFactor
import blang.core.LogScaleFactor
import java.util.List
import java.util.ArrayList
import static blang.runtime.internals.objectgraph.StaticUtils.node

class SplitSampler implements Sampler {
  @SampledVariable public Split split
  @ConnectedFactor public List<LogScaleFactor> numericFactors
  List<List<LogScaleFactor>> sortedNumericFactors
  
  override execute(Random rand) {
    throw new UnsupportedOperationException("TODO: auto-generated method stub")
  }
  
  def List<Double> normalizedTipInclusionProbabilities() {
    val result = new ArrayList
    for (index : 0 ..< split.nCells) {
      val p = eval(index, true)
      val q = eval(index, false)
      result.add(p / (p + q))
    }
    return result
  }
  
  def private double eval(int index, boolean value) {
    val indic = split.tipIndicators.get(index)
    val backup = indic.included
    indic.included = value
    val factors = sortedNumericFactors.get(index)
    var result = 0.0
    for (factor : factors)
      result += factor.logDensity
    indic.included = backup
    return result
  }
  
  override setup(SamplerBuilderContext context) {
    // TODO: check for loops
    sortedNumericFactors = new ArrayList
    for (indicator : split.tipIndicators) {
      val factors = context.connectedFactors(node(indicator))
      sortedNumericFactors.add(factors.toArray as List<LogScaleFactor>)
    }
    return true
  }
}
