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
import blang.mcmc.internals.ExponentiatedFactor
import org.eclipse.xtend.lib.annotations.Data

class SingleCellHMM implements TidilySerializable {
  
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
    
    val indices = data.positions.indices(chromosome) 
    var int len = indices.size 
    logReads = newDoubleArrayOfSize(len)
    logGCs = newDoubleArrayOfSize(len)
    for (Index<Integer> position : indices) {
      logReads.set(position.key, Math::log(data.readCounts.get(chromosome, position).intValue))
      logGCs.set(position.key, Math::log(data.gcContents.get(chromosome, position).doubleValue))
    }
    ChromoPostProcessor::addToPlot(logReads, logGCs) 
  }
  
  @Data
  static abstract class Base implements HMM {    
    protected val SingleCellHMM hmm
    protected val TransitionMatrix transition
    
    override transitionProbabilities(int t) {
      transition
    }
    
    override initialProbabilities() {
      blang.types.StaticUtils::fixedSimplex(ones(nStates) / nStates)
    }
    
    override length() {
      hmm.logReads.length
    } 
  }

  static class Standard extends Base { 
    new(SingleCellHMM hmm) {
      super(hmm, hmm.transitionMatrix)
    }
    override observationLogDensity(int t, int state) {
      hmm.readCountModel.logDensity(hmm.logGCs.get(t), hmm.logReads.get(t), state) 
    }
  }
  
  private def transitionMatrix() {
    val matrix = dense(nStates, nStates)
    val switchPr = switchProbability.doubleValue // TODO: change this to continuous time to have nice interpol
    if (switchPr < 0.0 || switchPr > 1.0) blang.types.StaticUtils::invalidParameter
    val offDiagonal = switchPr / (nStates - 1.0)
    matrix.editInPlace[x, y, v| if (x === y) 1.0 - switchProbability.doubleValue else offDiagonal]
    return blang.types.StaticUtils::fixedTransitionMatrix(matrix)
  }
    
  def double logMarginal() {
    val hmm = new Standard(this)
    val result = HMMComputations::logMarginalProbability(hmm, anneal)
    if (!anneal.present) return result
    return if (result == Double.NEGATIVE_INFINITY)
      ExponentiatedFactor::annealedMinusInfinity(anneal.get.doubleValue)
    else
      result
  }
  
  override serialize(Context context) {
    val hmm = new Standard(this)
    val samples = HMMComputations::sample(random, hmm)
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