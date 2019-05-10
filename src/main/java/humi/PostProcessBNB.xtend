package humi

import blang.inits.experiments.Experiment
import java.io.File
import blang.inits.Arg
import briefj.BriefIO
import blang.distributions.BetaNegativeBinomial
import bayonet.distributions.Random
import blang.types.StaticUtils
import java.util.List
import java.util.ArrayList

class PostProcessBNB extends Experiment {
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
    val rs = read("r.csv")
    val als = read("alpha.csv") 
    val bes = read("beta.csv") 
    for (i : 0 ..< rs.size) {
      dataset++
      
      // create a dataset
      var collected = 0
      val dist = BetaNegativeBinomial::distribution(
        StaticUtils::fixedReal(rs.get(i)),
        StaticUtils::fixedReal(als.get(i)),
        StaticUtils::fixedReal(bes.get(i))
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