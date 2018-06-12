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

@Samplers(CorruptGibbsSampler)
class CorruptPhylo {
  @Accessors(PUBLIC_GETTER)
  val PerfectPhylo reconstruction 
  val CellLocusMatrix tipInclPrs
  
  // Initialize with star tree
  new (CellLocusMatrix tipInclPrs) {
    this(new PerfectPhylo(tipInclPrs.cells, tipInclPrs.loci), tipInclPrs)
    if (!fixedTips)
      throw new RuntimeException("Caching might not work.")
  }
  
  new (PerfectPhylo phylo, CellLocusMatrix tipInclPrs) {
    if (phylo.cells != tipInclPrs.cells || 
        phylo.loci != tipInclPrs.loci)
      throw new RuntimeException
    this.tipInclPrs = tipInclPrs
    this.reconstruction = phylo
  }
  
  // un-annealed
  var double _cache = Double.NaN
  def logProbability() {
    ensureCache 
    return _cache
  }
  
  private def void ensureCache() {
    if (Double.isNaN(_cache))
      resetCache
  }
  
  def void resetCache() {
    _cache = 0.0
    for (locus : loci) 
      _cache += logProbability(locus)
  }
  
  def logProbability(Locus locus) {
    var sum = 0.0
    val tips = reconstruction.getTips(locus)
    for (entry : tips.entrySet) {
      val included = entry.value
      val pr = tipInclPrs.getTipAsDouble(entry.key, locus)
      sum += if (included) Math.log(pr) else Math.log1p(- pr)
    }
    return sum
  }
  
  def boolean fixedTips() { return tipInclPrs instanceof ReadOnlyCLMatrix }
  
  def loci()  { tipInclPrs.loci }
  def cells() { tipInclPrs.cells}
  
  def void priorSample(Random rand) {
    reconstruction.sampleUniform(rand)
    resetCache
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
    resetCache
    for (locus : loci)
      _gibbSample(rand, annealingParameters, locus)
  }
  
  // private as doing only one locus retriggers full likelihood computation
  private def void _gibbSample(Random rand, double annealingParameter, Locus locus) {
    ensureCache
    val likelihoodBefore = logProbability(locus)
    reconstruction.tree.collapseEdge(locus) 
    SplitSampler::sampleInPlace(reconstruction.tree, locus, cellInclusionLogProbabilities(annealingParameter, locus), rand)
    val likelihoodAfter = logProbability(locus)
    _cache += - likelihoodBefore + likelihoodAfter
  }
  
  static val LOG_EPSILON = -1e6 // we rely on division to cancel things, so keep non-zero
  static val LOG_ONE_MINUS_EPSILON = Math.log1p(Math.exp(LOG_EPSILON))
  def Map<Cell,SubtreeLikelihood> cellInclusionLogProbabilities(double annealingParameter, Locus locus) {
    val result = new LinkedHashMap
    for (cell : cells) {
      if (annealingParameter == 0.0) {
        result.put(cell, SubtreeLikelihood::missingTip) 
      } else {
        val inclPr = tipInclPrs.getTipAsDouble(cell, locus)
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
}