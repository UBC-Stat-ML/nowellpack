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
import blang.inits.Arg
import blang.inits.DefaultValue
import blang.inits.DesignatedConstructor
import blang.inits.ConstructorArg
import blang.inits.Implementations

/**
 * Creates various multi-resolution HMMs for efficient custom tempering
 */
class SingleCellHMMs implements TidilySerializable {
  
  static class Configs {
    @Arg  @DefaultValue("8")
    public int nStates = 8
    
    @Arg                                @DefaultValue("1")
    public Random demarginalizationRandom = new Random(1)
    
    val Optional<AnnealingParameter> annealingParameter
    
    val AnnealingStrategy annealingStrategy
    
    @DesignatedConstructor
    new(@ConstructorArg("annealingStrategy") AnnealingStrategy annealingStrategy) {
      annealingParameter = if (annealingStrategy instanceof Exponentiation) Optional.empty else Optional.of(new AnnealingParameter)
      this.annealingStrategy = annealingStrategy
    }
  }
  
  @Implementations(Exponentiation, Prefix, MultiLevel)
  static interface AnnealingStrategy {
    def double logMarginal(SingleCellHMMs enclosing)
  }
  
  static class Exponentiation implements AnnealingStrategy {
    @Arg   @DefaultValue("1")
    public int thinning = 1
    override double logMarginal(SingleCellHMMs enclosing) {
      return HMMComputations::logMarginalProbability(enclosing.getHMM(thinning))
    }
  }
  
  static class Prefix implements AnnealingStrategy {
    override logMarginal(SingleCellHMMs enclosing) {
      return HMMComputations::logMarginalProbability(enclosing.getHMM(1), Optional.of(enclosing.configs.annealingParameter.get))
    }
  }
  
  static class MultiLevel implements AnnealingStrategy {
    @Arg @DefaultValue("4")
    public      int b = 4
    @Arg @DefaultValue("3")
    public      int n = 3
    override logMarginal(SingleCellHMMs enclosing) {
      val beta = enclosing.configs.annealingParameter.get.doubleValue
      if (beta == 0.0) return 0.0
      if (beta == 1.0) return HMMComputations::logMarginalProbability(enclosing.getHMM(1))
      val int l = Math.floor(beta * n) as int
      val int u = l+1
      val log_pi_l = if (l == 0) 0.0 else HMMComputations::logMarginalProbability(enclosing.getHMM(thinning(l)))
      val log_pi_u = HMMComputations::logMarginalProbability(enclosing.getHMM(thinning(u)))
      val double lambda = n * beta - l
      return (1.0 - lambda) * log_pi_l + lambda * log_pi_u
    }
    private def int thinning(int i) {
      if (i == 0) throw new RuntimeException
      return Math::pow(b, n - i) as int
    }
  }
  
  val Configs configs
    
  @SkipDependency(isMutable = false)
  val double [] logReads
  
  @SkipDependency(isMutable = false)
  val double [] logGCs
  
  val ReadCountModel readCountModel
  val RealVar switchRate
  
  new(SingleCellData data, Index<String> chromosome, ReadCountModel readCountModel, RealVar switchRate, Configs configs) {
    ChromoPostProcessor::nStates = configs.nStates // hack    
    this.configs = configs
    this.readCountModel = readCountModel
    this.switchRate = switchRate
    
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
  static class SingleCellHMM implements HMM {    
    val SingleCellHMMs hmm
    val TransitionMatrix transition
    val int thinning
    
    override transitionProbabilities(int t) {
      transition
    }
    
    override initialProbabilities() {
      blang.types.StaticUtils::fixedSimplex(ones(hmm.configs.nStates) / hmm.configs.nStates)
    }
    
    override length() {
      hmm.logGCs.size / thinning
    } 
    
    override observationLogDensity(int t, int state) {
      hmm.readCountModel.logDensity(hmm.logGCs.get(t * thinning), hmm.logReads.get(t * thinning), state) 
    }
  }
  
  def getHMM(int thinning) {
    new SingleCellHMM(this, jcTransitionMatrix(configs.nStates, deltaTimeMu(thinning)), thinning)
  }
  
  def double deltaTimeMu(int thinning) {
    val rate = switchRate.doubleValue
    if (rate < 0.0) blang.types.StaticUtils::invalidParameter
    return rate * thinning
  }

  static def jcTransitionMatrix(int size, double deltaTimeMu) {
    val matrix = dense(size, size)
    if (deltaTimeMu < 0) throw new RuntimeException
    val diag = 1.0/size + (size - 1.0)/size*Math::exp(-deltaTimeMu)
    val offd = 1.0/size - 1.0/size*Math::exp(-deltaTimeMu)
    matrix.editInPlace[x, y, v| if (x === y) diag else offd]
    return blang.types.StaticUtils::fixedTransitionMatrix(matrix)
  }
    
  def double logMarginal() {
    val result = configs.annealingStrategy.logMarginal(this)
    if (configs.annealingStrategy instanceof Exponentiation) return result
    return if (result == Double.NEGATIVE_INFINITY)
      ExponentiatedFactor::annealedMinusInfinity(configs.annealingParameter.get.doubleValue)
    else
      result
  }
  
  override serialize(Context context) {
    val samples = HMMComputations::sample(configs.demarginalizationRandom, getHMM(1))
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