package chromobreak

import hmm.HMM
import blang.core.RealVar
import blang.types.TransitionMatrix

import static xlinear.MatrixOperations.*
import blang.types.Index
import blang.runtime.internals.objectgraph.SkipDependency
import hmm.HMMComputations
import blang.inits.experiments.tabwriters.TidilySerializable
import blang.inits.experiments.tabwriters.TidySerializer.Context
import bayonet.distributions.Random
import blang.types.AnnealingParameter
import java.util.Optional
import blang.runtime.Runner
import java.nio.file.Files

class SingleCellHMM implements HMM, TidilySerializable {
  
  public static val int nStates = 8
  
  val Optional<AnnealingParameter> anneal
    
  @SkipDependency(isMutable = false)
  val double [] logReads
  
  @SkipDependency(isMutable = false)
  val double [] logGCs
  
  @SkipDependency(isMutable = false)
  val Random random = new Random(1)
  
  val ReadCountModel readCountModel
  val RealVar switchProbability
  
  new(SingleCellData data, Index<String> chromosome, ReadCountModel readCountModel, RealVar switchProbability, Optional<AnnealingParameter> anneal) {
    this.anneal = anneal
    this.readCountModel = readCountModel
    this.switchProbability = switchProbability
    
    var int len = data.positions.indices(chromosome).size 
    logReads = newDoubleArrayOfSize(len)
    logGCs = newDoubleArrayOfSize(len)
    for (Index<Integer> position : data.positions.indices(chromosome)) {
      logReads.set(position.key, Math::log(data.readCounts.get(chromosome, position).intValue))
      logGCs.set(position.key, Math::log(data.gcContents.get(chromosome, position).doubleValue))
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
    val switchPr = switchProbability.doubleValue // TODO: change this to continuous time to have nice interpol
    if (switchPr < 0.0 || switchPr > 1.0) blang.types.StaticUtils::invalidParameter
    val offDiagonal = switchPr / (nStates - 1.0)
    matrix.editInPlace[x, y, v| if (x === y) 1.0 - switchProbability.doubleValue else offDiagonal]
    transition = blang.types.StaticUtils::fixedTransitionMatrix(matrix)
  }
    
  def double logMarginal() {
    preprocess()
    return HMMComputations::logMarginalProbability(this, anneal)
  }
  
  override serialize(Context context) {
    val samples = HMMComputations::sample(random, this)
    for (t : 0 ..< samples.size)
      context.recurse(samples.get(t), "positions", t)
  }
  
  def static void main(String [] args) {
    val tempDirWithPrefix = Files.createTempDirectory("SingleCellHMM")
    val runner = Runner::create(tempDirWithPrefix.toFile,
      "--model", "chromobreak.SingleCell",
      "--model.data.source", "data/SA922_uncor_gc_simplified.csv",
      "--model.data.gcContents.name", "value",
      "--model.data.readCounts.name", "value" ,
      "--model.data.readCounts.dataSource", "data/SA922_uncor_reads_split/5.csv",
      "--model.useDataAnneal", "true"
    )
    val model = runner.model as SingleCell
    model.components
    for (hmm : model.hmms.values)
      println(hmm.logMarginal)
  }
}