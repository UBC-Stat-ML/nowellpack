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
import corrupt.post.ConcatenationMatrix
import java.util.Set
import java.util.Collection

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
    val cache = pickCacheImpl(this, tipInclPrs) 
    this.cache = if (corrupt.CorruptPhylo.testCacheCorrectness) 
      new DebugCache(cache, new NoCache(this))
    else
      cache
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
    return logProbability(locus, tips, tipInclPrs)
  }
  
  def static logProbability(Locus locus, Map<Cell, Boolean> tips, CellLocusMatrix tipInclPrs) {
    var sum = 0.0
    for (entry : tips.entrySet) {
      val included = entry.value
      val pr = tipInclPrs.get(entry.key, locus)
      if (pr == Double::NEGATIVE_INFINITY) return Double::NEGATIVE_INFINITY // convention: if BinaryCLMatrix gets into invalid state, return - INF; then just bail out
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
  
  def void sample(Random rand, double annealingParameter) {
    sample(rand, annealingParameter, shuffledLoci(rand), cells)
  }
  
  def void sampleWithoutCellReallocation(Random rand, double annealingParameter) {
    sample(rand, annealingParameter, shuffledLoci(rand), Collections::emptyList)
  }
  
  def void sample(Random rand, double annealingParameter, Collection<Locus> sampledLoci, Collection<Cell> sampledCells) {
    for (locus : sampledLoci) {
      reconstruction.tree.collapseEdge(locus) 
      SplitSampler::sampleInPlace(reconstruction.tree, locus, inclusionLogProbabilities(annealingParameter, locus), rand)
    }
    for (cell : sampledCells) {
      reconstruction.tree.collapseEdge(cell)
      CellSampler::sampleInPlace(reconstruction.tree, cell, inclusionLogProbabilities(annealingParameter, cell), rand)
    }
    cache.reset
  }

  def Map<Cell,SubtreeLikelihood> inclusionLogProbabilities(double annealingParameter, Locus locus) {
    return inclusionLogProbabilities(annealingParameter, locus, cells, tipInclPrs)
  }
  def Map<Locus,SubtreeLikelihood> inclusionLogProbabilities(double annealingParameter, Cell cell) {
    return inclusionLogProbabilities(annealingParameter, cell, loci, tipInclPrs)
  }
  def static <T> Map<T,SubtreeLikelihood> inclusionLogProbabilities(
    double annealingParameter, 
    TreeNode reference, 
    Collection<T> orthogonals,
    CellLocusMatrix tipInclPrs
  ) {
    val result = new LinkedHashMap<T,SubtreeLikelihood>
    for (orthogonal : orthogonals) {
      if (annealingParameter == 0.0) {
        result.put(orthogonal, SubtreeLikelihood::missingTip) 
      } else {
        val inclPr = switch reference {
          Locus : tipInclPrs.get(orthogonal as Cell, reference)
          Cell :  tipInclPrs.get(reference, orthogonal as Locus)
          default : throw new RuntimeException
        }
        var logP = annealingParameter * Math.log(inclPr)
        var logQ = annealingParameter * Math.log1p(- inclPr)
             if (logP < LOG_EPSILON) { logP = LOG_EPSILON; logQ = LOG_ONE_MINUS_EPSILON }
        else if (logQ < LOG_EPSILON) { logQ = LOG_EPSILON; logP = LOG_ONE_MINUS_EPSILON }
        result.put(orthogonal, SubtreeLikelihood::tip(logP, logQ))
      }
    }
    return result
  }
  static val LOG_EPSILON = -1e6 // we rely on division to cancel things, so keep non-zero
  static val LOG_ONE_MINUS_EPSILON = Math.log1p(Math.exp(LOG_EPSILON))
  
  override toString(){ reconstruction.toString }
  
  // Caching stuff
    
  static interface Cache {
    def Set<Locus> loci()
    def double cachedLogPr()
    def boolean initialized()
    def void reset()
  }
  
  def private static Cache pickCacheImpl(CorruptPhylo phylo, CellLocusMatrix matrix) {
    return switch (matrix) {
      ConcatenationMatrix : new ConcatenationCache(phylo, matrix)
      ReadOnlyCLMatrix : new ReadOnlyMatrixCache(phylo, matrix)
      NoisyBinaryCLMatrix : new NoisyBinaryCache(phylo, matrix)
      default : {
        BriefLog::warnOnce("No cache found for " + matrix.class.simpleName + "... this may make certain infer schemes very slow")
        new NoCache(phylo) 
      }
    }
  }
  
  public static boolean testCacheCorrectness = false
  
  private static class ConcatenationCache implements Cache {
    val List<Cache> caches
    new (CorruptPhylo phylo, ConcatenationMatrix concatenation) {
      caches = new ArrayList
      for (matrix : concatenation.matrices) 
        caches.add(pickCacheImpl(phylo, matrix))
    }
    override loci() { throw new RuntimeException }
    override cachedLogPr() {
      var sum = 0.0
      for (cache : caches) 
        sum += cache.cachedLogPr
      return sum
    }
    override initialized() {
      for (cache : caches)
        if (!cache.initialized)
          return false
      return true
    }
    override reset() {
      for (cache : caches) cache.reset
    }
  }
  
  private static class ReadOnlyMatrixCache implements Cache {
    val CorruptPhylo phylo
    val ReadOnlyCLMatrix matrix
    new(CorruptPhylo phylo, ReadOnlyCLMatrix matrix) { 
      this.phylo = phylo
      this.matrix = matrix
    }
    var double cachedValue = Double.NaN
    override loci() { matrix.loci }
    override cachedLogPr() { cachedValue }
    override initialized() { !Double.isNaN(cachedValue) }
    override reset() {
      cachedValue = 0.0
      for (locus : matrix.loci) {
         val tips = phylo.reconstruction.getTips(locus)
         cachedValue += logProbability(locus, tips, matrix)
      }
    }
  }
  
  private static class NoisyBinaryCache implements Cache {
    val CorruptPhylo phylo
    val NoisyBinaryCLMatrix noisyMatrix
    var List<NoiseStatistics> stats = null
    new(CorruptPhylo phylo, NoisyBinaryCLMatrix matrix) { 
      this.phylo = phylo
      this.noisyMatrix = matrix
    }
    override loci() { noisyMatrix.loci }
    override cachedLogPr() {
      noisyMatrix.sumLogPrs(stats) 
    }
    override initialized() {
      return stats !== null
    }
    override reset() {
      stats = new ArrayList
      for (locus : noisyMatrix.loci) {
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
  }
  
  private static class NoCache implements Cache {
    val CorruptPhylo phylo
    new(CorruptPhylo phylo) { this.phylo = phylo }
    override loci() { phylo.loci }
    override cachedLogPr() {
      var sum = 0.0
      for (locus : phylo.loci)
        sum += phylo.logProbability(locus)
      return sum
    }
    override initialized() { true }
    override reset() {}
  }
  
  private static class DebugCache implements Cache {
    val Cache c1
    val Cache c2
    new(Cache c1, Cache c2) { 
      this.c1 = c1
      this.c2 = c2
    }
    override loci() { throw new RuntimeException }
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
  }
}