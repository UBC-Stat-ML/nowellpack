package humi

import blang.types.Plated
import blang.types.Index
import blang.types.StaticUtils
import blang.core.RealVar


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
}