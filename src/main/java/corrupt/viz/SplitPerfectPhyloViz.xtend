package corrupt.viz

import blang.inits.Arg
import blang.inits.experiments.Experiment
import corrupt.PerfectPhylo
import corrupt.post.ReadOnlyCLMatrix
import corrupt.viz.PerfectPhyloViz
import java.util.List
import viz.core.PublicSize
import java.util.Optional
import corrupt.GenomeMap

/**
 * Split by chromosomes for large matrices
 */
class SplitPerfectPhyloViz extends Experiment {
  
  @Arg public PerfectPhylo phylo 
  @Arg public List<ReadOnlyCLMatrix> matrices
  @Arg PublicSize size
  @Arg Optional<PerfectPhylo> ref
  @Arg Optional<List<Integer>> colourCodes
  @Arg Optional<String> suffix
  
  override run() {
    val output = results.child("output")
    val allLoci = PerfectPhyloViz::allLoci(matrices)
    val map = new GenomeMap(allLoci)
    for (chr : map.orderedChromosomes) {
      val viz = new PerfectPhyloViz(phylo, matrices, size, ref, colourCodes, Optional.of(map.orderedLoci(chr).toSet))
      val outFile = output.getFileInResultFolder("chr_" + chr + (if (suffix.present) "_" + suffix else "") + ".pdf") 
      viz.output(outFile)
    }
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}