package corrupt

import blang.inits.experiments.Experiment
import blang.inits.Arg
import java.io.File
import static corrupt.post.CLMatrixUtils.*


import blang.inits.DefaultValue


class Distances extends Experiment {
  
  @Arg File reference
  @Arg File guess
  
  
  @Arg @DefaultValue("null.csv")
  public File SNVref
  
  @Arg @DefaultValue("null.csv")
  public File SNVprob
  
  @Arg @DefaultValue("null.csv")
  public File Evalref
  
  
  @Arg  @DefaultValue("0.5")
  public double threshold = 0.5
  
  
  override run() {
    val refPhylo   = PerfectPhylo::parseNewick(reference)
    val guessPhylo = PerfectPhylo::parseNewick(guess)

    
	if (SNVref.exists() && SNVprob.exists() && Evalref.exists()) {
      val refSNV = fromCSV(SNVref, true)
      val guessSNV = fromCSV(SNVprob, true)
      val evalRef = fromCSV(Evalref, true)

      results.getTabularWriter("distances") => [
        write("metric" -> "distance", "value" -> distance(refPhylo, guessPhylo))
        write("metric" -> "locusEdgeDistance", "value" -> locusEdgeDistance(refPhylo, guessPhylo))
        write("metric" -> "sensitivity", "value" -> SNVSensitivity(round(refSNV), round(guessSNV), evalRef)) 
        write("metric" -> "specificity", "value" -> SNVSpecificity(round(refSNV), round(guessSNV), evalRef))
        write("metric" -> "Youden-index", "value" -> SNVSensitivity(round(refSNV), round(guessSNV), evalRef) + SNVSpecificity(round(refSNV), round(guessSNV), evalRef) -1)
        write("metric" -> "F1-score", "value" -> SNVF1(round(refSNV), round(guessSNV), evalRef))  
      ]  
    } else{
    	  results.getTabularWriter("distances") => [
        write("metric" -> "distance",          "value" -> distance(refPhylo, guessPhylo))
        write("metric" -> "locusEdgeDistance", "value" -> locusEdgeDistance(refPhylo, guessPhylo))
        
      ]
    }
  }
  
  public static def void main(String [] args) {
    Experiment.start(args)
  }
}
