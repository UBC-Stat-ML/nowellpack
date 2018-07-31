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
  @Arg public int nCells
  @Arg public int nLoci
  
  @Arg              @DefaultValue("1")
  val Random treeRand = new Random(1)
  @Arg              @DefaultValue("1")
  val Random dataRand = new Random(1)
  
  @Arg(description = "Controls how hard the problem is; lower value are easier to infer")
          @DefaultValue("0.3")
  public double stdDev = 0.3
  
  @Arg           @DefaultValue("false")
  public boolean useFPFNRates = false
  
  @Arg    @DefaultValue("0.01")
  public double fpRate = 0.01
  
  @Arg    @DefaultValue("0.05")
  public double fnRate = 0.05
  
  public static val TREE_FILE = "phylo.newick"
  public static val DATA_FILE = "tipInclusionProbabilities.csv"
  
  override run() {
    // generate and write phylogen
    val phylo = new PerfectPhylo(syntheticCells(nCells), syntheticLoci(nLoci))
    phylo.sampleUniform(treeRand)
    BriefIO::write(results.getFileInResultFolder(TREE_FILE), phylo.toNewick)
  
    // generate tips
    val data = 
      if (useFPFNRates)
        CLMatrixUtils::syntheticPerturbedBinaryMatrix(dataRand, phylo, fpRate, fnRate)
      else
        CLMatrixUtils::syntheticInclusionPrs(dataRand, phylo, stdDev)
    data.toCSV(results.getFileInResultFolder(DATA_FILE), CLMatrixUtils::fromPhylo(phylo)) 
    
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