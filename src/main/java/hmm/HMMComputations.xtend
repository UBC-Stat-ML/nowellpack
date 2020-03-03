package hmm

import static extension xlinear.MatrixExtensions.*
import static xlinear.MatrixOperations.*
import xlinear.Matrix
import java.util.List
import java.util.ArrayList
import java.util.Optional
import blang.types.AnnealingParameter

class HMMComputations {
  
  def static double logMarginalProbability(HMM hmm) {
    return logMarginalProbability(hmm, Optional.empty)
  }
  
  /**
   * Log probability of observations summing over latent states
   */
  def static double logMarginalProbability(HMM hmm, Optional<AnnealingParameter> anneal) {
    return backward(hmm, null, anneal)
  }
  
  /**
   * Or null if no path have positive probability
   */
  def static List<Integer> sample(bayonet.distributions.Random random, HMM hmm) {
    val savedBackwardVectors = new ArrayList<Matrix>
    if (backward(hmm, savedBackwardVectors, Optional.empty) == Double.NEGATIVE_INFINITY) return null
    val result = new ArrayList<Integer>
    for (t : 0 ..< savedBackwardVectors.size) {
      val currentVector = savedBackwardVectors.get(savedBackwardVectors.size-t-1)
      if (t === 0) {
        val init = hmm.initialProbabilities
        for (s : 0 ..< currentVector.nEntries)
          currentVector.set(s, currentVector.get(s) * init.get(s))
      } else {
        val transition = hmm.transitionProbabilities(t - 1)
        val prevState = result.get(t - 1)
        for (s : 0 ..< currentVector.nEntries)
          currentVector.set(s, currentVector.get(s) * transition.get(prevState, s))
      }
      currentVector /= currentVector.sum
      val sample = random.nextCategorical(currentVector.vectorToArray)
      result.add(sample)
    }
    return result
  }
  
  private def static double backward(HMM hmm, List<Matrix> savedBackwardVectors, Optional<AnnealingParameter> anneal) {
    val double beta = if (anneal.present) anneal.get.doubleValue else 1.0
    val boolean saveVectors = savedBackwardVectors !== null
    
    val lengthToConsider = Math.ceil(beta * hmm.length) as int
    val partialAnneal = lengthToConsider - beta * hmm.length
    val indexToAnneal = if (anneal.present) (lengthToConsider - 1) else -1
    if (lengthToConsider == 0) return 0.0
    var initial = hmm.initialProbabilities
    val lastStepSize = if (lengthToConsider == 1) initial.nEntries else hmm.transitionProbabilities(lengthToConsider - 2).nCols
    var vector = ones(lastStepSize)
    var sumLogs = 0.0
    for (t : (lengthToConsider - 1) .. 0) {
      var foundPossibleState = false
      for (s : 0 ..< vector.nEntries) {
          var obsLogDensity = hmm.observationLogDensity(t, s) 
          
          if (Double.isNaN(obsLogDensity) || obsLogDensity == Double.POSITIVE_INFINITY) 
            throw new RuntimeException
            
          if (obsLogDensity !== Double.NEGATIVE_INFINITY)
            foundPossibleState = true
          
          if (indexToAnneal == t) {
            if (partialAnneal == 0.0)
              /*
               * this is needed in (rare but occurring) case where 
               * partialAnneal == 0.0 and obsLogDensity = -INF
               */ 
              obsLogDensity = 0 
            else
              obsLogDensity *= partialAnneal
          }
          
          vector.set(s, Math::log(vector.get(s)) + obsLogDensity)
        }
        
      if (!foundPossibleState) {
        if (saveVectors) 
          throw new RuntimeException
        else 
          return Double.NEGATIVE_INFINITY
      }
            
      val logNorm = expNormalize(vector)
      
      if (saveVectors) savedBackwardVectors.add(vector.copy)
      if (logNorm == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY
      sumLogs += logNorm
      
      if (t > 0)
        vector = hmm.transitionProbabilities(t - 1) * vector

    }
    return  sumLogs + Math.log(vector.dot(initial))
  }
  
  def static double expNormalize(Matrix logProbs) {
    var max = Double.NEGATIVE_INFINITY
    for(i : 0 ..< logProbs.nEntries)
      max = Math.max(max, logProbs.get(i))
    for(i : 0 ..< logProbs.nEntries)
      logProbs.set(i,Math.exp(logProbs.get(i)-max))
    val sum = logProbs.sum
    logProbs /= sum
    return max + Math.log(sum)
  }
}