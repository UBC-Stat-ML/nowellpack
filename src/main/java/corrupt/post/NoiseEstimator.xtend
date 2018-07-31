package corrupt.post

import corrupt.PerfectPhylo
import java.io.File
import blang.inits.experiments.Experiment
import blang.inits.Arg

/**
 * Inputs:
 *  - zero-one data matrix
 *  - tree samples
 * 
 * Output: estimates of false positive, negative rates.
 * 
 * Note similarity with ScoreTree, but here we want want one of the objects to 
 * be loadable from a data matrix, also the role of what is iterable/truth is 
 * changed, and finally we want a more precise notion of distance to distinguish 
 * FP and FN. 
 */
class NoiseEstimator extends Experiment {
  @Arg public File data
  @Arg public File tree
  
  /**
   * Assume phylo represents the "truth" and binaryMatrix is a noisy observed version of 
   * the induced latent clade indicators. What are the false positive and 
   * false negative rates?
   */
  static def void estimate(PerfectPhylo phylo, CellLocusMatrix binaryMatrix, NoiseStatistics statistics) {
    if (phylo.cells != binaryMatrix.cells || phylo.loci != binaryMatrix.loci)
      throw new RuntimeException
    for (locus : phylo.loci) 
      statistics.add(locus, phylo.getTips(locus), binaryMatrix)
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
  
  override run() {
    val loadedData = CLMatrixUtils::fromCSV(data)
    val loadedTree = PerfectPhylo::parseNewick(tree)
    val noise = new NoiseStatistics
    estimate(loadedTree, loadedData, noise)
    println("FP = " + noise.fpRate)
    println("FN = " + noise.fnRate)
    results.getTabularWriter("errorRates") => [
      write("type" -> "FP", "value" -> noise.fpRate)
      write("type" -> "FN", "value" -> noise.fnRate)
    ]
  }
}