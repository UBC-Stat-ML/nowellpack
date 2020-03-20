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
import blang.core.IntVar
import hmm.SparseTransitionMatrix

/**
 * Creates various multi-resolution HMMs for efficient custom tempering
 */
class SingleCellHMMs implements TidilySerializable {
  
  static class Configs {
    @Arg          @DefaultValue("true")
    public boolean bufferState = true
    
    @Arg                                @DefaultValue("1")
    public Random demarginalizationRandom = new Random(1)
    
    @Arg(description = "Can be set to infinity but useful to set finite max to control resource usage")          
            @DefaultValue("20")
    public int maxStates = 20
    
    val Optional<AnnealingParameter> annealingParameter
    
    val AnnealingStrategy annealingStrategy
    
    @DesignatedConstructor
    new(@ConstructorArg("annealingStrategy") @DefaultValue("Exponentiation") AnnealingStrategy annealingStrategy) {
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
    @Arg @DefaultValue("10")
    public      int b = 10
    @Arg @DefaultValue("2")
    public      int n = 2
    @Arg            @DefaultValue("1")
    public      int baseThinning = 1
    @DesignatedConstructor
    new () {
      MultiLevelPT::configs = this
    }
    override logMarginal(SingleCellHMMs enclosing) {
      val beta = enclosing.configs.annealingParameter.get.doubleValue
      if (beta == 0.0) return 0.0
      if (beta == 1.0) return HMMComputations::logMarginalProbability(enclosing.getHMM(baseThinning))
      val int l = Math.floor(beta * n) as int
      val int u = l+1
      val log_pi_l = if (l == 0) 0.0 else HMMComputations::logMarginalProbability(enclosing.getHMM(baseThinning*thinning(l)))
      val log_pi_u = HMMComputations::logMarginalProbability(enclosing.getHMM(baseThinning*thinning(u)))
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
  
  val IntVar _nStates
  
  def int nStates() {
    val result = _nStates.intValue
    if (result < 2 || result > configs.maxStates) blang.types.StaticUtils::invalidParameter
    return result + if (configs.bufferState) 1 else 0
  }
  
  val ReadCountModel readCountModel
  val RealVar switchRate
  
  new(SingleCellData data, Index<String> chromosome, ReadCountModel readCountModel, RealVar switchRate, Configs configs, IntVar nStates) {
    this.configs = configs
    this._nStates = nStates
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
      blang.types.StaticUtils::fixedSimplex(ones(hmm.nStates) / hmm.nStates)
    }
    
    override length() {
      hmm.logGCs.size / thinning
    } 
    
    override observationLogDensity(int t, int state) {
      if (hmm.configs.bufferState) {
        if (state === hmm.nStates - 1)
          return -Math.log(10)
      }
      hmm.readCountModel.logDensity(hmm.logGCs.get(t * thinning), hmm.logReads.get(t * thinning), state) 
    }
  }
  
  def getHMM(int thinning) {
    val deltaTimeMu = deltaTimeMu(thinning)
    return new SingleCellHMM(
      this, 
      (if (configs.bufferState) bufferedTransitionMatrix(nStates, deltaTimeMu) else jcTransitionMatrix(nStates, deltaTimeMu)), 
      thinning
    )
  }
  
  def double deltaTimeMu(int thinning) {
    val rate = switchRate.doubleValue
    if (rate < 0.0) blang.types.StaticUtils::invalidParameter
    return rate * thinning
  }
  
  static def bufferedTransitionMatrix(int size, double deltaTimeMu) {
    val matrix = sparse(size, size)
    if (deltaTimeMu < 0) throw new RuntimeException
    val prNoChange = Math::exp(-deltaTimeMu)
    val bufferState = size - 1
    for (s : 0 .. bufferState-1) {
      matrix.set(s, s, prNoChange)
      matrix.set(s, bufferState, 1.0 - prNoChange)
      matrix.set(bufferState, s, 1.0 / (size-1))
    }
    return new SparseTransitionMatrix(matrix)
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
    val bufferState = if (configs.bufferState) nStates-1 else -1
    for (t : 0 ..< samples.size)
      context.recurse(processState(samples.get(t), bufferState), "positions", t)
  }
  
  def static processState(int value, int bufferState) {
    if (value === bufferState) "NA"
    else "" + value
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