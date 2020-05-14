package corrupt


import blang.inits.experiments.Experiment
import blang.inits.Arg
import corrupt.post.CellLocusMatrix
import blang.inits.DefaultValue
import briefj.BriefIO
import bayonet.distributions.Random
import java.util.Collections
import java.util.ArrayList




import static extension briefj.BriefIO.write

class GrowTree extends Experiment {
  @Arg public CellLocusMatrix matrix
  @Arg public PerfectPhylo phylo
 
  @Arg        @DefaultValue("1")
  public Random rand = new Random(1);
  
  @Arg @DefaultValue("2")
  int   nIterations = 2
  
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
  	
    var out = BriefIO.output(results.getFileInResultFolder("SNV_prob.csv"))
    
    val shuffledLoci = new ArrayList(matrix.loci)
    Collections::shuffle(shuffledLoci, rand)
    out.append("cells," + "loci," + "snv_prob" + "\n") 
    for (i : 0 ..< nIterations)
      for (locus : shuffledLoci) {
        if (phylo.tree.nodes.contains(locus)) {
          phylo.tree.collapseEdge(locus)
        }
        val likelihoods = CorruptPhylo::inclusionLogProbabilities(1.0, locus, matrix.cells, matrix)
        if (i==0){
        		SplitSampler::maximizeInPlace(phylo.tree, locus, likelihoods, false, null)
        	
        	}
        	else{
        		SplitSampler::maximizeInPlace(phylo.tree, locus, likelihoods, true, out)	
        	}
      }
    out.close()  
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
