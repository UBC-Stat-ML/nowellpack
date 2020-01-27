package chromobreak

import briefj.collections.Counter
import org.eclipse.xtend.lib.annotations.Data
import static java.lang.Math.*
import bayonet.math.NumericalUtils
import static blang.types.StaticUtils.*
import java.util.List

@Data
class PoissonGamma {
  
  val Counter<MultiplierVector> vectors = new Counter
  val int size
  val double shape
  val double rate
  
  def void withVector(MultiplierVector v) {
    if (!vectors.empty) throw new RuntimeException
    vectors.incrementCount(v, 1.0)
  }
  
  def double logProbability(List<Integer> observations) {
    if (observations.size != size) throw new RuntimeException
    val n = observations.stream.mapToInt[intValue].sum
    val logNorm = log(vectors.totalCount)
    val alpha = shape
    val beta = rate
    var logFactorials = 0.0
    for (observation : observations)
      logFactorials += logFactorial(observation)
    var logSum = NEGATIVE_INFINITY
    for (vector : vectors) {
      var current = log(vectors.getCount(vector)) - logNorm
      var sum = 0.0
      for (i : 0 ..< size) {
        val mu = vector.get(i)
        sum += mu
        current += observations.get(i) * log(mu)
      }
      current -= (n + alpha) * log(beta + sum)
      logSum = NumericalUtils::logAdd(logSum, current)
    }
    return logSum + alpha * log(beta) + logGamma(n + alpha) - logGamma(alpha) - logFactorials
  }
  
  static interface MultiplierVector {
    def double get(int i)
  }
}