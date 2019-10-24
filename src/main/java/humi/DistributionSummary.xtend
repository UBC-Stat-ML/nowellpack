package humi

import blang.core.IntDistribution
import blang.types.Index
import humi.CountFrequencies
import briefj.collections.Counter
import bayonet.math.NumericalUtils
import humi.Monitor
import blang.types.Plated
import java.util.function.Supplier
import blang.core.IntVar
import blang.core.RealVar
import blang.distributions.BetaNegativeBinomial
import blang.types.StaticUtils
import blang.distributions.YuleSimon
import blang.distributions.NegativeBinomialMeanParam
import blang.distributions.Poisson
import bayonet.distributions.Random
import blang.distributions.Generators

class DistributionSummary {
  
  static val cutoff = 20
  static val winsorizationP = 0.9 // 0.99 seems to lead to long tail which gets computationally costy
  static val guard = 10_000
  
  val static rand = new Random(1)
  
  
  def static registerMonitors(
    Plated<Monitor> visibleCloneNumbers, 
    Plated<Monitor> truncatedMeans,
    Plated<Monitor> truncatedSqMeans,
    Plated<Monitor> winsorizedMeans,
    Plated<Monitor> conditionalWinsorizedMeans,
    Index<Integer> target,
    Index<String> experiment,
    Supplier<IntDistribution> dist,
    Plated<IntVar> initialPopCounts,
    RealVar p1,
    RealVar p2
  ) {
    
    val initialPopCount = initialPopCounts.get(target)
    
    // register visible clone number monitors
    visibleCloneNumbers.get(target, experiment).init(
      if (p2 === null) {
        val RealVar lambda = p1
        new RealVar() {
          override doubleValue() {
            val p0 = Math.exp(dist.get.logDensity(0))
            return (1.0 - p0) * initialPopCount.intValue * lambda.doubleValue
          }
        }
      } else {
        val RealVar shape = p1
        val RealVar rate = p2
        new RealVar() {
          override doubleValue() {
            val p0 = Math.exp(dist.get.logDensity(0))
            val lambda = Generators::gamma(rand, shape.doubleValue, rate.doubleValue) 
            return (1.0 - p0) * initialPopCount.intValue * lambda
          }
        }
      }
    )
    
    // register truncated mean monitors
    if (!truncatedMeans.get(target).initialized)
      truncatedMeans.get(target).init(new RealVar() {
        override doubleValue() {
          DistributionSummary::mean(DistributionSummary::truncatedNormalizedCounter(dist.get))
        }
      })
      
    if (!truncatedSqMeans.get(target).initialized)
      truncatedSqMeans.get(target).init(new RealVar() {
        override doubleValue() {
          DistributionSummary::meanSq(DistributionSummary::truncatedNormalizedCounter(dist.get))
        }
      })
    
    // winsorized mean
    if (!winsorizedMeans.get(target).initialized)
      winsorizedMeans.get(target).init(new RealVar() {
        override doubleValue() {
          return winsorizedMean(dist.get, winsorizationP)
        }
      })
    
    // conditional version
    if (!conditionalWinsorizedMeans.get(target).initialized) 
      conditionalWinsorizedMeans.get(target).init(new RealVar() {
        override doubleValue() {
          return conditionalWinsorizedMean(dist.get, winsorizationP)
        }
      })
  }
  
  def static double conditionalWinsorizedMean(IntDistribution distribution, double p) {
    if (p < 0.5 || p > 1.0) throw new RuntimeException
    var sum = 0.0
    var x = 0
    val normalization = 1.0 - Math::exp(distribution.logDensity(0))
    while (sum/normalization < p) {
      if (x > guard) return Double.NaN
      x++  
      sum += Math::exp(distribution.logDensity(x))
    }
    val cutOff = x
    var result = 0.0
    var mass = 0.0
    for (y : 0 ..< cutOff) {
      val currentMass = Math::exp(distribution.logDensity(y))
      mass += currentMass
      result += y * currentMass
    }
    result += cutOff * (1.0 - mass)
    return result
  }
  
  def static double winsorizedMean(IntDistribution distribution, double p) {
    val Counter<Integer> truncated = truncate(distribution, p)
    if (truncated === null) return Double.NaN
    return mean(truncated)
  }
  
  def private static Counter<Integer> truncate(IntDistribution distribution, double p) {
    val result = new Counter
    var mass = 0.0
    var c = 0
    while (true) {
      if (c > guard) return null
      
      val cur =  Math::exp(distribution.logDensity(c))
      
      if (mass + cur < p) {
        result.setCount(c, cur)
        mass += cur
      } else {
        result.setCount(c, 1.0 - mass)
        NumericalUtils::checkIsClose(result.totalCount, 1.0)
        return result
      }
      c++
    }
  }
  
  def static Counter<Integer> truncatedNormalizedCounter(CountFrequencies frequencies) {
    val result = new Counter
    for (c : 1 .. cutoff)
      if (frequencies.distinctCounts.contains(c))
        result.setCount(c, frequencies.frequency(c))
    result.normalize
    return result
  }
  
  def static Counter<Integer> truncatedNormalizedCounter(IntDistribution dist) {
    val result = new Counter
    for (c : 1 .. cutoff)
      result.setCount(c, Math.exp(dist.logDensity(c)))
    result.normalize
    return result
  }
  
  // input: a normalized int valued dist, stored in a counter
  def static double mean(Counter<Integer> counter, int degree) {
    NumericalUtils::checkIsClose(1.0, counter.totalCount) 
    var sum = 0.0
    for (c : counter.keySet)
      sum += (c ** degree) * counter.getCount(c)
    return sum
  }
  
  def static double mean(Counter<Integer> counter) { mean(counter, 1) }
  def static double meanSq(Counter<Integer> counter) { mean(counter, 2) }
  
  def static Supplier<IntDistribution> bnb(RealVar r, RealVar alpha, RealVar beta) {
    new Supplier<IntDistribution>() {
      override IntDistribution get() {
        return BetaNegativeBinomial::distribution(r, alpha, beta)
      }
    }
  }
  
  def static Supplier<IntDistribution> nb(RealVar mean, RealVar od) {
    new Supplier<IntDistribution>() {
      override IntDistribution get() {
        return NegativeBinomialMeanParam::distribution(mean, od)
      }
    }
  }
  
  def static Supplier<IntDistribution> ys(RealVar rho) {
    new Supplier<IntDistribution>() {
      override IntDistribution get() {
        return YuleSimon::distribution(rho)
      }
    }
  }
  
  def static Supplier<IntDistribution> poi(RealVar mean) {
    new Supplier<IntDistribution>() {
      override IntDistribution get() {
        return Poisson::distribution(mean)
      }
    }
  }
  
  def static Supplier<IntDistribution> mixBnb(RealVar pi, RealVar r1, RealVar alpha1, RealVar beta1, RealVar r2, RealVar alpha2, RealVar beta2) {
    new Supplier<IntDistribution>() {
      override IntDistribution get() {
        return IntMixture::distribution({
                if (pi.doubleValue < 0.0 || pi.doubleValue > 1.0) StaticUtils::invalidParameter
                StaticUtils::fixedSimplex(#[pi.doubleValue, 1.0 - pi.doubleValue])
              }, #[ 
                BetaNegativeBinomial::distribution(r1, alpha1, beta1), 
                BetaNegativeBinomial::distribution(r2, alpha2, beta2)
              ]
            )
      }
    }
  }
  
  def static Supplier<IntDistribution> mixNb(RealVar pi, RealVar mean1, RealVar od1, RealVar mean2, RealVar od2) {
    new Supplier<IntDistribution>() {
      override IntDistribution get() {
        return IntMixture::distribution({
                if (pi.doubleValue < 0.0 || pi.doubleValue > 1.0) StaticUtils::invalidParameter
                StaticUtils::fixedSimplex(#[pi.doubleValue, 1.0 - pi.doubleValue])
              }, #[ 
                NegativeBinomialMeanParam::distribution(mean1, od1), 
                NegativeBinomialMeanParam::distribution(mean2, od2)
              ]
            )
      }
    }
  }
  
  def static Supplier<IntDistribution> mixYs(RealVar pi, RealVar rho1, RealVar rho2) {
    new Supplier<IntDistribution>() {
      override IntDistribution get() {
        return IntMixture::distribution({
                if (pi.doubleValue < 0.0 || pi.doubleValue > 1.0) StaticUtils::invalidParameter
                StaticUtils::fixedSimplex(#[pi.doubleValue, 1.0 - pi.doubleValue])
              }, #[ 
                YuleSimon::distribution(rho1), 
                YuleSimon::distribution(rho2)
              ]
            )
      }
    }
  }
}