package humi

import blang.distributions.Poisson
import static blang.types.StaticUtils.*
import blang.core.IntDistribution
import static java.lang.Math.*
import blang.core.RealVar
import java.util.Random

class CensoredExchangeableUtils {
  
  static var maxIters = 1_000_000
  def static int sampleTotalCount(Random rand, IntDistribution pmf, RealVar poissonRate, int nObserved) {
    val uniform = rand.nextDouble
    val logQ = pmf.logDensity(0)
    val q = Math::exp(logQ)
    val log1MQ = Math.log1p(-q)
    val logNorm = Poisson::distribution(fixedReal((1.0 - q) * poissonRate.doubleValue)).logDensity(nObserved)
    val poi = Poisson::distribution(poissonRate)
    var cum = 0.0
    for (n : nObserved ..< nObserved + maxIters) {
      cum += exp(poi.logDensity(n) + logBinomial(n, nObserved) + (n - nObserved) * logQ + nObserved * log1MQ - logNorm)
      if (cum >= uniform) return n
    }
    throw new RuntimeException
  }
  
  def static void main(String [] args) {
    val q = 0.2
    val nObserved = 2
    val lam = 5.4
    val poi = Poisson::distribution(fixedReal(lam))
    
    
    println("norm = " + Math::exp(Poisson::distribution(fixedReal((1.0 - q) * lam)).logDensity(nObserved)))
    
    var sum = 0.0
    for (n : nObserved ..< 10) {
      sum += exp(poi.logDensity(n)) * exp(logBinomial(n, nObserved)) * pow(q, n - nObserved) * pow(1.0-q, nObserved)
      println(sum)
    }
  }
}