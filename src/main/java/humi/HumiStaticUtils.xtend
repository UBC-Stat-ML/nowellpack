package humi

import blang.types.Plated
import blang.types.Index
import blang.types.StaticUtils
import blang.core.RealVar
import java.util.List
import java.util.ArrayList
import blang.distributions.NegativeBinomialMeanParam
import blang.core.IntDistribution

class HumiStaticUtils { 
  
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
  
  def static nbMix(List<RealVar> means, List<RealVar> overds) {
    val result = new ArrayList<IntDistribution>
    for (i : 0 ..< means.size)
      result.add(NegativeBinomialMeanParam::distribution(means.get(i), overds.get(i)))
    return result
  }
}