package humi

import blang.types.Plated
import blang.types.Index
import blang.types.StaticUtils
import blang.core.RealVar
import java.util.List
import java.util.ArrayList
import blang.distributions.NegativeBinomialMeanParam
import blang.core.IntDistribution
import blang.core.IntVar

import static bayonet.math.SpecialFunctions.*

class HumiStaticUtils {
  
  def static List<IntVar> controlIndicatorsList(HumiData data, Plated<IntVar> controlIndicators) {
    // extract indicators for the controls
    val indicators = new ArrayList<IntVar>
    for (Index<String> gene : data.genes.indices.filter[data.isControl(it)]) {
      for (Index<Integer> target : data.targets.indices(gene)) {
        indicators.add(controlIndicators.get(target))
      }
    }
    return indicators
  }
  
  def static double logPoissonCensoringFormula(double poissonRate, int numberGreaterThanZero, double probabilityForOneZero) {
    return 
        numberGreaterThanZero * Math::log(poissonRate) 
      + (probabilityForOneZero - 1.0) * poissonRate 
      - StaticUtils::logFactorial(numberGreaterThanZero)
  }
  
  def static RealVar orOneIfControl(Plated<RealVar> targets, String controlPrefix, Index<String> gene, Index<Integer> target) {
    val geneName = gene.key
    if (geneName.startsWith(controlPrefix))
      return StaticUtils::fixedReal(1.0)
    else
      return targets.get(gene, target)
  }
  
  def static nbMix(RealVar rate, List<RealVar> means, List<RealVar> overds) {
    val result = new ArrayList<IntDistribution>
    for (i : 0 ..< means.size)
      result.add(NegativeBinomialMeanParam::distribution([means.get(i).doubleValue * rate.doubleValue], overds.get(i)))
    return result
  }
}