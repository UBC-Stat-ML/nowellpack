package corrupt

import briefj.Indexer
import corrupt.post.CellLocusMatrix
import java.util.Map
import java.util.LinkedHashMap
import bayonet.distributions.Random
import blang.mcmc.Samplers
import org.eclipse.xtend.lib.annotations.Accessors

@Samplers(CorruptGibbsSampler)
class CorruptPhylo {
  @Accessors(PUBLIC_GETTER)
  val PerfectPhylo phylo 
  val CellLocusMatrix tipInclPrs
  
  new (CellLocusMatrix tipInclPrs) {
    this.tipInclPrs = tipInclPrs
    this.phylo = new PerfectPhylo(tipInclPrs.cells, tipInclPrs.loci)
  }

  // un-annealed
  def logProbability() {
    var sum = 0.0
    for (locus : loci) {
      val tips = phylo.getTips(locus)
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
    phylo.sampleUniform(rand)
  }

  def void gibbSample(Random rand, double annealingParameter) {
    for (locus : loci) 
      _gibbSample(rand, annealingParameter, locus)
  }
  
  // private as doing only one locus retriggers full likelihood computation
  private def void _gibbSample(Random rand, double annealingParameter, Locus locus) {
    phylo.tree.collapseEdge(locus) 
    SplitSampler::sampleInPlace(phylo.tree, locus, cellInclusionLogProbabilities(annealingParameter, locus), rand)
  }
  
  def Map<Cell,SubtreeLikelihood> cellInclusionLogProbabilities(double annealingParameter, Locus locus) {
    val result = new LinkedHashMap
    for (cell : cells) {
      val inclPr = tipInclPrs.getTipAsDouble(cell, locus)
      val logP = annealingParameter * Math.log(inclPr)
      val logQ = annealingParameter * Math.log(1.0 - inclPr)
      result.put(cell, SubtreeLikelihood::tip(logP, logQ))
    }
    return result
  }
  
  override toString(){ phylo.toString }
}