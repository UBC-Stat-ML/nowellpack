package hmm

import blang.types.TransitionMatrix
import blang.types.Simplex

interface HMM {
  /**
   * from t to t+1
   */
  def TransitionMatrix transitionProbabilities(int t)
  
  def Simplex initialProbabilities()
  
  def int length()
  
  def double observationLogDensity(int t, int state)
}