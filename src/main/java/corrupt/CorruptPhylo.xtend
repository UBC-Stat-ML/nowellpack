package corrupt

import corrupt.post.CellLocusMatrix
import java.util.Map
import java.util.LinkedHashMap
import bayonet.distributions.Random
import blang.mcmc.Samplers
import org.eclipse.xtend.lib.annotations.Accessors
import java.util.ArrayList
import java.util.Collections
import corrupt.post.ReadOnlyCLMatrix
import java.util.List
import bayonet.math.NumericalUtils
import briefj.BriefLog
import corrupt.post.NoiseStatistics
import corrupt.post.NoisyBinaryCLMatrix
import blang.runtime.internals.objectgraph.SkipDependency

@Samplers(CorruptGibbsSampler)
class CorruptPhylo {
  @Accessors(PUBLIC_GETTER)
  @SkipDependency(isMutable = true) // declared in the conditioning of LogPot; otherwise wrong behaviour
  val PerfectPhylo reconstruction 
  
  @SkipDependency(isMutable = true)
  val CellLocusMatrix tipInclPrs
  
  @Accessors(PUBLIC_GETTER)
  @SkipDependency(isMutable = true)
  val Cache cache
  
  // Initialize with star tree
  new (CellLocusMatrix tipInclPrs) {
    this(new PerfectPhylo(tipInclPrs.cells, tipInclPrs.loci), tipInclPrs)
  }
  
  new (PerfectPhylo phylo, CellLocusMatrix tipInclPrs) {
    if (phylo.cells != tipInclPrs.cells || 
        phylo.loci != tipInclPrs.loci)
      throw new RuntimeException
    this.tipInclPrs = tipInclPrs
    this.reconstruction = phylo
    this.cache = pickCacheImpl(this, tipInclPrs.class) 
  }
  
  // un-annealed
  def logProbability() {
    ensureCache 
    return cache.cachedLogPr + tipInclPrs.logNormalization 
  }
  
  private def void ensureCache() {
    if (!cache.initialized)
      cache.reset
  }
  
  def logProbability(Locus locus, Map<Cell, Boolean> tips) {
    var sum = 0.0
    for (entry : tips.entrySet) {
      val included = entry.value
      val pr = tipInclPrs.get(entry.key, locus)
      if (pr == Double::NEGATIVE_INFINITY) return Double::NEGATIVE_INFINITY
      sum += if (included) Math.log(pr) else Math.log1p(- pr)
    }
    return sum
  }
  
  def logProbability(Locus locus) {
    val tips = reconstruction.getTips(locus)
    return logProbability(locus, tips)
  }
  
  def boolean fixedTips() { return tipInclPrs instanceof ReadOnlyCLMatrix }
  
  def loci()  { tipInclPrs.loci }
  def cells() { tipInclPrs.cells}
  
  def void priorSample(Random rand) {
    reconstruction.sampleUniform(rand)
    cache.reset
  }
  
  private def List<Locus> shuffledLoci(Random rand) {
    val shuffled = new ArrayList(loci)
    Collections::shuffle(shuffled, rand) 
    return shuffled
  }
  
  var List<Locus> _shuffled = null
  var int _index = 0
  def void nextGibbs(Random rand, double annealingParameter) {
    val index = _index % loci.size
    if (index === 0)
      _shuffled = shuffledLoci(rand)
    _gibbSample(rand, annealingParameter, _shuffled.get(index))
    _index++
  }
  
  def void gibbsTest(Random rand, double annealingParameters) {    
    cache.reset
    for (locus : loci)
      _gibbSample(rand, annealingParameters, locus)
  }
  
  // private as doing only one locus retriggers full likelihood computation
  private def void _gibbSample(Random rand, double annealingParameter, Locus locus) {
    ensureCache
    val tipsBefore = reconstruction.getTips(locus)
    reconstruction.tree.collapseEdge(locus) 
    SplitSampler::sampleInPlace(reconstruction.tree, locus, cellInclusionLogProbabilities(annealingParameter, locus), rand)
    val tipsAfter = reconstruction.getTips(locus)
    cache.update(locus, tipsBefore, tipsAfter)
  }
  
  static val LOG_EPSILON = -1e6 // we rely on division to cancel things, so keep non-zero
  static val LOG_ONE_MINUS_EPSILON = Math.log1p(Math.exp(LOG_EPSILON))
  def Map<Cell,SubtreeLikelihood> cellInclusionLogProbabilities(double annealingParameter, Locus locus) {
    val result = new LinkedHashMap
    for (cell : cells) {
      if (annealingParameter == 0.0) {
        result.put(cell, SubtreeLikelihood::missingTip) 
      } else {
        val inclPr = tipInclPrs.get(cell, locus)
        var logP = annealingParameter * Math.log(inclPr)
        var logQ = annealingParameter * Math.log1p(- inclPr)
             if (logP < LOG_EPSILON) { logP = LOG_EPSILON; logQ = LOG_ONE_MINUS_EPSILON }
        else if (logQ < LOG_EPSILON) { logQ = LOG_EPSILON; logP = LOG_ONE_MINUS_EPSILON }
        result.put(cell, SubtreeLikelihood::tip(logP, logQ))
      }
    }
    return result
  }
  
  override toString(){ reconstruction.toString }
  
  // Caching stuff
    
  static interface Cache {
    def double cachedLogPr()
    def boolean initialized()
    def void reset()
    def void update(Locus locus, Map<Cell, Boolean> tipsBefore, Map<Cell, Boolean> tipsAfter)
  }
  
  def private static Cache pickCacheImpl(CorruptPhylo phylo, Class<? extends CellLocusMatrix> matrixType) {
    val Cache selection = switch (matrixType) {
      case ReadOnlyCLMatrix : new ReadOnlyMatrixCache(phylo)
      case NoisyBinaryCLMatrix : new NoisyBinaryCache(phylo)
      default : {
        BriefLog::warnOnce("No cache found for " + matrixType.simpleName + "... this may make certain infer schemes very slow")
        new NoCache(phylo)
      }
    }
    if (corrupt.CorruptPhylo.testCacheCorrectness) 
      return new DebugCache(selection, new NoCache(phylo))
    else
      return selection
  }
  
  public static boolean testCacheCorrectness = false
  
  private static class ReadOnlyMatrixCache implements Cache {
    val CorruptPhylo phylo
    new(CorruptPhylo phylo) { this.phylo = phylo }
    var double cachedValue = Double.NaN
    override cachedLogPr() { cachedValue }
    override initialized() { !Double.isNaN(cachedValue) }
    override reset() {
      cachedValue = 0.0
      for (locus : phylo.loci) 
        cachedValue += phylo.logProbability(locus)
    }
    override update(Locus locus, Map<Cell, Boolean> tipsBefore, Map<Cell, Boolean> tipsAfter) {
      val loglBefore = phylo.logProbability(locus, tipsBefore)
      val loglAfter = phylo.logProbability(locus, tipsAfter)
      cachedValue += - loglBefore + loglAfter
    }
  }
  
  private static class NoisyBinaryCache implements Cache {
    val CorruptPhylo phylo
    val NoisyBinaryCLMatrix noisyMatrix
    var List<NoiseStatistics> stats = null
    new(CorruptPhylo phylo) { 
      this.phylo = phylo
      this.noisyMatrix = phylo.tipInclPrs as NoisyBinaryCLMatrix
    }
    override cachedLogPr() {
      noisyMatrix.sumLogPrs(stats) 
    }
    override initialized() {
      return stats !== null
    }
    override reset() {
      stats = new ArrayList
      for (locus : phylo.loci) {
        val currentStat = 
          if (noisyMatrix.global) {
            if (stats.empty) {
              val global = new NoiseStatistics
              stats.add(global)
            }
            stats.get(0)
          } else {
            val specific = new NoiseStatistics
            stats.add(specific)
            specific
          }
        currentStat.add(locus, phylo.reconstruction.getTips(locus), noisyMatrix.binaryMatrix) 
      }
    }
    override update(Locus locus, Map<Cell, Boolean> tipsBefore, Map<Cell, Boolean> tipsAfter) {
      stats.get(noisyMatrix.parameter(locus)).subtract(locus, tipsBefore, noisyMatrix.binaryMatrix) 
      stats.get(noisyMatrix.parameter(locus)).add(locus, tipsAfter, noisyMatrix.binaryMatrix) 
    }
  }
  
  private static class NoCache implements Cache {
    val CorruptPhylo phylo
    new(CorruptPhylo phylo) { this.phylo = phylo }
    override cachedLogPr() {
      var sum = 0.0
      for (locus : phylo.loci)
        sum += phylo.logProbability(locus)
      return sum
    }
    override initialized() { true }
    override reset() {}
    override update(Locus locus, Map<Cell, Boolean> tipsBefore, Map<Cell, Boolean> tipsAfter) {}
  }
  
  private static class DebugCache implements Cache {
    val Cache c1
    val Cache c2
    new(Cache c1, Cache c2) { 
      this.c1 = c1
      this.c2 = c2
    }
    override cachedLogPr() {
      val v1 = c1.cachedLogPr
      val v2 = c2.cachedLogPr
      if (Double.isInfinite(v1) && Double.isInfinite(v2))
        return v1
      NumericalUtils::checkIsClose(v1, v2, 1e-10)
      return v1
    }
    override initialized() {
      return c1.initialized && c2.initialized
    }
    override reset() {
      c1.reset
      c2.reset
    }
    override update(Locus locus, Map<Cell, Boolean> tipsBefore, Map<Cell, Boolean> tipsAfter) {
      c1.update(locus, tipsBefore, tipsAfter)
      c2.update(locus, tipsBefore, tipsAfter)
    }
  }
}