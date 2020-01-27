package hmm

import static extension xlinear.MatrixExtensions.*
import static xlinear.MatrixOperations.*
import xlinear.Matrix

class HMMComputation {
  
  def static double logMarginalProbability(HMM hmm) {
    val len = hmm.length
    if (len == 0) return 0.0
    var initial = hmm.initialProbabilities
    val lastStepSize = if (len == 1) initial.nEntries else hmm.transitionProbabilities(len - 2).nCols
    var vector = ones(lastStepSize)
    var sumLogs = 0.0
    for (t : (len - 1) .. 0) {
      for (s : 0 ..< vector.nEntries) {
        val obsLogDensity = hmm.observationLogDensity(t, s) 
        
        if (Double.isNaN(obsLogDensity)) 
          throw new RuntimeException
          
        vector.set(s, Math::log(vector.get(s)) + obsLogDensity)
      }
      
      val logNorm = expNormalize(vector)
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