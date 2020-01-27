package chromobreak

import hmm.HMM
import blang.core.RealVar
import blang.types.TransitionMatrix
import hmm.HMMComputation

import static xlinear.MatrixOperations.*
import blang.types.Index
import blang.runtime.internals.objectgraph.SkipDependency

class SingleCellHMM implements HMM {
  
  val int nStates = 6
  
  val int thinning = 10
    
  @SkipDependency(isMutable = false)
  val double [] logReads
  
  @SkipDependency(isMutable = false)
  val double [] logGCs
  
  val ReadCountModel readCountModel
  val RealVar switchProbability
  
  new(SingleCellData data, Index<String> chromosome, ReadCountModel readCountModel, RealVar switchProbability) {
    this.readCountModel = readCountModel
    this.switchProbability = switchProbability
    
    var int len = 0 
    for (Index<Integer> position : data.positions.indices(chromosome).filter[key % thinning == 0]) {
      len++
    }
    
    logReads = newDoubleArrayOfSize(len)
    logGCs = newDoubleArrayOfSize(len)
    for (Index<Integer> position : data.positions.indices(chromosome).filter[key % thinning == 0]) {
      logReads.set(position.key / thinning, Math::log(data.readCounts.get(chromosome, position).intValue))
      logGCs.set(position.key / thinning, Math::log(data.gcContents.get(chromosome, position).doubleValue))
    }
    ChromoPostProcessor::addToPlot(logReads, logGCs) 
  }
  
  TransitionMatrix transition = null
  
  override transitionProbabilities(int t) {
    transition
  }
  
  override initialProbabilities() {
    blang.types.StaticUtils::fixedSimplex(ones(nStates) / nStates)
  }
  
  override length() {
    logReads.length
  }
  
  override observationLogDensity(int t, int state) {
    readCountModel.logDensity(logGCs.get(t), logReads.get(t), state) 
  }
  
  private def preprocess() {
    val matrix = dense(nStates, nStates)
    val switchPr = switchProbability.doubleValue
    if (switchPr < 0.0 || switchPr > 1.0) blang.types.StaticUtils::invalidParameter
    val offDiagonal = switchPr / (nStates - 1.0)
    matrix.editInPlace[x, y, v| if (x === y) 1.0 - switchProbability.doubleValue else offDiagonal]
    transition = blang.types.StaticUtils::fixedTransitionMatrix(matrix)
  }
    
  def double logMarginal() {
    preprocess()
    return HMMComputation::logMarginalProbability(this)
  }
  
}