package corrupt

import blang.inits.experiments.Experiment
import blang.inits.Arg
import java.util.Collections



import static corrupt.CorruptStaticUtils.*
import blang.inits.DefaultValue
import bayonet.distributions.Random
import briefj.BriefIO
import corrupt.post.CLMatrixUtils



import static extension corrupt.post.CLMatrixUtils.toCSV
import static extension corrupt.post.CLMatrixUtils.toCSV_pt
import static extension corrupt.post.CLMatrixUtils.toCSV_snv


class GenerateData extends Experiment {
  @Arg public int nCells
  @Arg public int nCNVLoci
  @Arg public int nSNVLoci

  
  
  @Arg              @DefaultValue("1")
  val Random treeRand = new Random(1)
  
  @Arg    @DefaultValue("1")
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
  
  
  @Arg @DefaultValue("0.01")
  public double seqError = 0.01
  
  @Arg @DefaultValue("0.01")
  public double seqDO = 0.01
  
  @Arg @DefaultValue("0.07")
  public double coverage = 0.07
  
  @Arg @DefaultValue("0.2")
  public double snvLociPercentage = 0.2
  
  @Arg @DefaultValue("5")
  public int cMax = 5
  
  
  
  public static val TREE_FILE = "phylo.newick"
  public static val TREE_FILE_1 = "phylo_1.newick"
  public static val MUT_REF_FILE = "SNV_ref.csv"
  public static val EVAL_REF_FILE = "SNV_eval_ref.csv"
  public static val DATA_FILE = "tipInclusionProbabilities.csv"
  public static val DATA_FILE_0 = "tipInclusionProbabilities_0.csv" //cnv tip
  public static val DATA_FILE_1 = "tipInclusionProbabilities_1.csv" //subset of snv tip
  public static val DATA_FILE_2 = "tipInclusionProbabilities_2.csv" // all of snv_tip
  
  
  override run() {
    // generate and write phylogen
    //val phylo = new PerfectPhylo(syntheticCells(nCells), syntheticLoci(nLoci, 0))
    val phylo = new PerfectPhylo(syntheticCells(nCells), syntheticLoci(nCNVLoci+ nSNVLoci))
    var shuffledLociId = newDoubleArrayOfSize(nCNVLoci+ nSNVLoci)
    for (locus:phylo.loci){
    	 shuffledLociId.set(locus.getIntegerId(),locus.getIntegerId())
    }
    Collections::shuffle(shuffledLociId, treeRand)
    for (locus:phylo.loci){
    	  var index = shuffledLociId.get(locus.getIntegerId()) 
    	  if (index < nCNVLoci){
    	  	locus.setType('cnv')
    	  	locus.setPrintType('0') 
    	  } else {
    	  	locus.setType('snv') //round or floor
    	  	if (index < (nCNVLoci + Math.round(snvLociPercentage*nSNVLoci))){ 
    	  	  locus.setPrintType('1') 
    	  	} else{
    	  	  locus.setPrintType('2')
    	  	}	  	
    	  }
    }  
    phylo.sampleUniform(treeRand)
    
    val referenceMutData = CLMatrixUtils::SyntheticReferenceMatrix(dataRand, phylo)
    referenceMutData.toCSV_snv(results.getFileInResultFolder(MUT_REF_FILE), true) 
    // generate tips
    val data = 
      if (useFPFNRates)
        CLMatrixUtils::syntheticPerturbedBinaryMatrix(dataRand, phylo, fpRate, fnRate, seqError, seqDO, coverage, cMax, results.getFileInResultFolder(EVAL_REF_FILE))
      else
        CLMatrixUtils::syntheticInclusionPrs(dataRand, phylo, stdDev, seqError, seqDO, coverage, cMax, results.getFileInResultFolder(EVAL_REF_FILE))
    
    data.toCSV(results.getFileInResultFolder(DATA_FILE), CLMatrixUtils::fromPhylo(phylo)) 

    data.toCSV_pt(results.getFileInResultFolder(DATA_FILE_0), '0', CLMatrixUtils::fromPhylo(phylo)) 
    data.toCSV_pt(results.getFileInResultFolder(DATA_FILE_1), '1', CLMatrixUtils::fromPhylo(phylo)) 
    data.toCSV_pt(results.getFileInResultFolder(DATA_FILE_2), '2', CLMatrixUtils::fromPhylo(phylo)) 
    

    // print probabilities
    val corruptPhylo = new CorruptPhylo(phylo, data)
    val likelihood = corruptPhylo.logProbability
    val prior = - CorruptStaticUtils::logNPerfectPhylo(nCells, nCNVLoci+ nSNVLoci)
    results.getTabularWriter("prs") => [
      write("type" -> "likelihood", "value" -> likelihood)
      write("type" -> "prior", "value" -> prior)
      write("type" -> "joint", "value" -> prior + likelihood)
    ]

    BriefIO::write(results.getFileInResultFolder(TREE_FILE), phylo.toNewick)
    
    var phylo_1 = phylo.collapseSubset("2")
    BriefIO::write(results.getFileInResultFolder(TREE_FILE_1), phylo_1.toNewick) 
  }
  
  static def void main(String [] args) {
    Experiment.startAutoExit(args)
  }
}