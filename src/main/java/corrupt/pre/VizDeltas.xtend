package corrupt.pre

import blang.inits.Arg
import blang.inits.experiments.Experiment
import blang.inits.experiments.tabwriters.factories.CSV
import corrupt.PerfectPhylo
import corrupt.post.CLMatrixUtils
import corrupt.post.ReadOnlyCLMatrix
import corrupt.viz.PerfectPhyloViz
import java.io.File
import java.util.ArrayList
import viz.core.Viz

class VizDeltas extends Experiment {
  
  @Arg public File tree
  
  @Arg public File matricesDir
  
  override run() {
    val phylo = PerfectPhylo::parseNewick(tree)
    val matrices= new ArrayList<ReadOnlyCLMatrix>
    for (var int i = 0; i < 3; i++) {
      val current = CSV::csvFile(matricesDir, "matrix-" + i)
      if (current !== null)
        matrices.add(ReadOnlyCLMatrix::readOnly(CLMatrixUtils::fromCSV(current)))
    }
    if (matrices.empty)
      throw new RuntimeException("No matrices found in " + matricesDir)
    PerfectPhyloViz::visualizePerChromosome(results.getFileInResultFolder("output"), phylo, matrices, Viz::fixHeight(300))
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}