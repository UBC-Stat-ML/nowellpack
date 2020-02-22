package chromobreak

import org.junit.Test
import hmm.HMM
import org.junit.Assert

import static extension xlinear.MatrixExtensions.*
import static xlinear.MatrixOperations.*
import blang.distributions.Normal
import com.google.common.base.Stopwatch
import java.util.concurrent.TimeUnit
import hmm.HMMComputations
import blang.validation.ExactInvarianceTest
import hmm.SimpleHMM

class TestHMMComputation {
  
  @Test
  def void testSampler() {
    val test = new ExactInvarianceTest()
    
    val builder = (new SimpleHMM.Builder)
      .setDynamics(blang.types.StaticUtils::fixedTransitionMatrix(
      #[
        #[0.3, 0.7],
        #[0.8, 0.2]
      ]
      ))
      .setEmissions(blang.types.StaticUtils::fixedTransitionMatrix(
      #[
        #[0.6, 0.4],
        #[0.8, 0.2]
      ]
      ))
      .setLength(4)
      .setInitial(blang.types.StaticUtils::fixedSimplex(0.4, 0.6))
    test.add(builder.build, [latents.get(0).intValue as double], [latents.get(1).intValue as double])
    test.check
  }
  
  @Test
  def void testLogPr() {
    val t0 = blang.types.StaticUtils::fixedTransitionMatrix(
      #[
        #[0.2, 0.7, 0.1],
        #[0.8, 0.1, 0.1]
      ]
    )
    val t1 = blang.types.StaticUtils::fixedTransitionMatrix(
      #[
        #[0.3, 0.7],
        #[0.8, 0.2],
        #[0.1, 0.9]
      ]
    )
    val e0 = blang.types.StaticUtils::fixedTransitionMatrix(
      #[
        #[0.3, 0.7],
        #[0.8, 0.2]
      ]
    )
    val e1 = blang.types.StaticUtils::fixedTransitionMatrix(
      #[
        #[0.3, 0.7],
        #[0.8, 0.2],
        #[0.1, 0.9]
      ]
    )
    val e2 = e0
    val i = blang.types.StaticUtils::fixedSimplex(0.4, 0.6)
    var sum = 0.0
    for (o0 : 0 .. 1) 
      for (o1 : 0 .. 1)
        for (o2 : 0 .. 1) {
          val hmm = new HMM() {
            override transitionProbabilities(int t) {
              return #[t0, t1].get(t)
            }
            override initialProbabilities() {
              return i
            }
            override length() {
              3
            }
            override observationLogDensity(int t, int state) {
              return Math::log(#[e0, e1, e2].get(t).get(state, #[o0, o1, o2].get(t)))
            }
          }
          sum += Math.exp(hmm.HMMComputations::logMarginalProbability(hmm))
        }
     Assert.assertEquals(1.0, sum, 1e-10)
  }
  
  @Test
  def void testOneObs() {
    
    val e0 = blang.types.StaticUtils::fixedTransitionMatrix(
      #[
        #[0.3, 0.7],
        #[0.8, 0.2]
      ]
    )
    
    val i = blang.types.StaticUtils::fixedSimplex(0.4, 0.6)
    var sum = 0.0

    for (o : 0 .. 1) {
      val hmm = new HMM() {
        override transitionProbabilities(int t) {
          throw new RuntimeException
        }
        override initialProbabilities() {
          return i
        }
        override length() {
          1
        }
        override observationLogDensity(int t, int state) {
          return Math::log(e0.get(state, o))
        }
      }
      sum += println(Math.exp(HMMComputations::logMarginalProbability(hmm)))
    }
     Assert.assertEquals(1.0, sum, 1e-10)
  }
  
  
  @Test
  def void performance() {
    val nStates = 5
    val transition = denseCopy(identity(nStates))
    val tr = blang.types.StaticUtils::fixedTransitionMatrix(transition)
    val len = 6000
    val gaussian = Normal::distribution(blang.types.StaticUtils::fixedReal(0.0), blang.types.StaticUtils::fixedReal(1.0))
    val hmm = new HMM() {
      
      override transitionProbabilities(int t) {
        return tr
      }
      
      override initialProbabilities() {
        blang.types.StaticUtils::fixedSimplex(ones(nStates) / nStates)
      }
      
      override length() {
        len
      }
      
      override observationLogDensity(int t, int state) {
        (gaussian.logDensity(state))
      }
      
    }
    println("started")
    HMMComputations::logMarginalProbability(hmm)
    val watch = Stopwatch.createStarted
    for (i : 0 .. 100)
      HMMComputations::logMarginalProbability(hmm)
    val elapsed = watch.elapsed(TimeUnit.MILLISECONDS)
    println("done in " + elapsed)
    Assert.assertTrue(elapsed < 6000)
  }
}