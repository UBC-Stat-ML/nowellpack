package humi

import blang.inits.experiments.Experiment
import java.io.File
import blang.inits.Arg
import briefj.BriefIO
import blang.distributions.YuleSimon
import bayonet.distributions.Random
import blang.types.StaticUtils
import java.util.List
import java.util.ArrayList

class PostProcessMix extends Experiment {
  @Arg File samples
  @Arg int thin 
  @Arg int burn 
  @Arg CountFrequencies freqs
  
  override run() {
    
    val output = results.getTabularWriter("samples")
    
    // check number of obs
    var nObserved = 0
    var dataset = 0
    for (count : freqs.distinctCounts) {
      val freq = freqs.frequency(count)
      nObserved += freq
      for (i : 0 ..< freq) 
        output.write("dataset" -> dataset, "count" -> count)
    }
    
    val rand = new Random(1)
    val pis = read("pi.csv")
    val rho1s = read("rho1.csv") 
    val rho2s = read("rho2.csv") 
    for (i : 0 ..< pis.size) {
      dataset++
      
      // create a dataset
      var collected = 0
      val dist = IntMixture::distribution(
            StaticUtils::fixedSimplex(pis.get(i), 1.0 - pis.get(i)), 
            #[
              YuleSimon::distribution(StaticUtils::fixedReal(rho1s.get(i))), 
              YuleSimon::distribution(StaticUtils::fixedReal(rho2s.get(i)))
            ]
          )
      while (collected < nObserved) {
        val sample = dist.sample(rand)
        if (sample > 0) {
          output.write("dataset" -> dataset, "count" -> sample)
          collected++
        }
      }
    }
  }
  
  def List<Double> read(String name) {
    val result = new ArrayList
    val file = new File(samples, name)
    for (line : BriefIO::readLines(file).indexCSV) {
      val iter = Integer.parseInt(line.get("sample"))
      if (iter >= burn && iter % thin === 0) {
        val v = Double.parseDouble(line.get("value"))
        result.add(v)
      }
    }
    return result
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
  
}