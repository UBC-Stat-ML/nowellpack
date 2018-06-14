package corrupt

import blang.inits.experiments.Experiment
import blang.inits.Arg

import static corrupt.CorruptStaticUtils.*
import blang.inits.DefaultValue
import bayonet.distributions.Random
import briefj.BriefIO
import corrupt.post.CLMatrixUtils

import static extension corrupt.post.CLMatrixUtils.toCSV

class GenerateData extends Experiment {
  @Arg int nCells
  @Arg int nLoci
  
  @Arg              @DefaultValue("1")
  val Random treeRand = new Random(1)
  @Arg              @DefaultValue("1")
  val Random dataRand = new Random(1)
  
  @Arg(description = "Controls how hard the problem is; lower value are easier to infer")
       @DefaultValue("0.3")
  double stdDev = 0.3
  
  override run() {
    // generate and write phylogen
    val phylo = new PerfectPhylo(syntheticCells(nCells), syntheticLoci(nLoci))
    phylo.sampleUniform(treeRand)
    BriefIO::write(results.getFileInResultFolder("phylo.newick"), phylo.toNewick)
  
    // generate tips
    val data = CLMatrixUtils::syntheticInclusionPrs(dataRand, phylo, stdDev)
    data.toCSV(results.getFileInResultFolder("tipInclusionProbabilities.csv"), CLMatrixUtils::fromPhylo(phylo)) 
    
    // print probabilities
    val corruptPhylo = new CorruptPhylo(phylo, data)
    val likelihood = corruptPhylo.logProbability
    val prior = - CorruptStaticUtils::logNPerfectPhylo(nCells, nLoci)
    results.getTabularWriter("prs") => [
      write("type" -> "likelihood", "value" -> likelihood)
      write("type" -> "prior", "value" -> prior)
      write("type" -> "joint", "value" -> prior + likelihood)
    ]
  }
  
  static def void main(String [] args) {
    Experiment.startAutoExit(args)
  }
}