package humi

import blang.inits.experiments.Experiment
import java.io.File
import blang.inits.Arg
import briefj.BriefIO
import blang.distributions.YuleSimon
import bayonet.distributions.Random
import blang.types.StaticUtils

class PostProcessSingle extends Experiment {
  @Arg File samples
  @Arg int thin 
  @Arg int burn 
  @Arg CountFrequencies freqs
  
  override run() {
    
    val output = results.getTabularWriter("samples")
    
    // check number of obs
    var nObserved = 0
    var dataset = 0
    for (count : freqs.counts) {
      val freq = freqs.frequency(count)
      nObserved += freq
      for (i : 0 ..< freq) 
        output.write("dataset" -> dataset, "count" -> count)
    }
    
    val rand = new Random(1)
    val rhos = new File(samples, "rho.csv")
    for (line : BriefIO::readLines(rhos).indexCSV) {
      val iter = Integer.parseInt(line.get("sample"))
      if (iter >= burn && iter % thin === 0) {
        dataset++
        val rho = Double.parseDouble(line.get("value"))
        // create a dataset
        var collected = 0
        val dist = YuleSimon::distribution(StaticUtils::fixedReal(rho))
        while (collected < nObserved) {
          val sample = dist.sample(rand)
          if (sample > 0) {
            output.write("dataset" -> dataset, "count" -> sample)
            collected++
          }
        }
      }
    }
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
  
}