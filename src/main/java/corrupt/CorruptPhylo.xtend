package corrupt

import corrupt.post.CellLocusMatrix
import java.util.Map
import java.util.LinkedHashMap
import bayonet.distributions.Random
import blang.mcmc.Samplers
import org.eclipse.xtend.lib.annotations.Accessors
import java.util.ArrayList
import java.util.Collections

@Samplers(CorruptGibbsSampler)
class CorruptPhylo {
  @Accessors(PUBLIC_GETTER)
  val PerfectPhylo reconstruction 
  val CellLocusMatrix tipInclPrs
  
  // TODO: cache the tip's logs if they are fixed
  // TODO: cache the likelihood based on tipInclPrs' hashcode?
  
  new (CellLocusMatrix tipInclPrs) {
    this.tipInclPrs = tipInclPrs
    this.reconstruction = new PerfectPhylo(tipInclPrs.cells, tipInclPrs.loci)
  }

  // un-annealed
  def logProbability() {
    var sum = 0.0
    for (locus : loci) {
      val tips = reconstruction.getTips(locus)
      for (entry : tips.entrySet) {
        val included = entry.value
        val pr = tipInclPrs.getTipAsDouble(entry.key, locus)
        sum += if (included) Math.log(pr) else Math.log(1.0 - pr)
      }
    }
    return sum
  }
  
  def loci()  { tipInclPrs.loci }
  def cells() { tipInclPrs.cells}
  
  def void priorSample(Random rand) {
    reconstruction.sampleUniform(rand)
  }

  def void gibbSample(Random rand, double annealingParameter) {
    val shuffled = new ArrayList(loci)
    Collections::shuffle(shuffled, rand) 
    for (locus : shuffled) 
      _gibbSample(rand, annealingParameter, locus)
  }
  
  // private as doing only one locus retriggers full likelihood computation
  private def void _gibbSample(Random rand, double annealingParameter, Locus locus) {
    reconstruction.tree.collapseEdge(locus) 
    SplitSampler::sampleInPlace(reconstruction.tree, locus, cellInclusionLogProbabilities(annealingParameter, locus), rand)
  }
  
  static val LOG_EPSILON = -1e6 // we rely on division to cancel things, so keep non-zero
  static val LOG_ONE_MINUS_EPSILON = Math.log1p(Math.exp(LOG_EPSILON))
  def Map<Cell,SubtreeLikelihood> cellInclusionLogProbabilities(double annealingParameter, Locus locus) {
    val result = new LinkedHashMap
    for (cell : cells) {
      val inclPr = tipInclPrs.getTipAsDouble(cell, locus)
      var logP = annealingParameter * Math.log(inclPr)
      var logQ = annealingParameter * Math.log1p(- inclPr)
           if (logP < LOG_EPSILON) { logP = LOG_EPSILON; logQ = LOG_ONE_MINUS_EPSILON }
      else if (logQ < LOG_EPSILON) { logQ = LOG_EPSILON; logP = LOG_ONE_MINUS_EPSILON }
      result.put(cell, SubtreeLikelihood::tip(logP, logQ))
    }
    return result
  }
  
  override toString(){ reconstruction.toString }
}