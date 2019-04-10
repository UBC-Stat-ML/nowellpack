package corrupt

import blang.inits.experiments.Experiment
import blang.inits.Arg
import corrupt.post.CellLocusMatrix
import blang.inits.DefaultValue

import static extension briefj.BriefIO.write
import java.util.Collections
import bayonet.distributions.Random
import java.util.ArrayList

class GrowTree extends Experiment {
  @Arg public CellLocusMatrix matrix
  @Arg public PerfectPhylo phylo
  
  @Arg @DefaultValue("2")
  int   nIterations = 2
  
  @Arg @DefaultValue("true")
  boolean randomize = true
  
  @Arg @DefaultValue("1")
  Random random = new Random(1)
  
  override run() {
    if (phylo.cells == matrix.cells)
      addLoci
    else if (phylo.loci == matrix.loci)
      addCells
    else
      throw new RuntimeException("The provided matrix and tree should either have (1) the same set of loci " + 
        "(in which cases extra cells are placed at the leaves) or (2) the same set of cells " + 
        "(in which case extra loci are placed in the tree)")
  
    results.getFileInResultFolder("grown.newick").write(phylo.toNewick)
  }
  
  def void addLoci() {
    for (i : 0 ..< nIterations) {
      val loci = new ArrayList(matrix.loci)
      if (randomize) 
        Collections.shuffle(loci, random)
      for (locus : loci) {
        if (phylo.tree.nodes.contains(locus)) {
          if (i == 0) throw new RuntimeException("Locus " + locus + " already placed in tree")
          phylo.tree.collapseEdge(locus)
        }
        val likelihoods = CorruptPhylo::inclusionLogProbabilities(1.0, locus, matrix.cells, matrix)
        SplitSampler::maximizeInPlace(phylo.tree, locus, likelihoods)
      }
    }
  }
  
  def void addCells() {
    for (cell : matrix.cells) {
      if (phylo.tree.nodes.contains(cell)) { 
        throw new RuntimeException("Cell " + cell + " already placed in tree")
      }
      val likelihoods = CorruptPhylo::inclusionLogProbabilities(1.0, cell, matrix.loci, matrix)
      CellSampler::maximizeInPlace(phylo.tree, cell, likelihoods)
    }
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}