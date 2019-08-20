package humi.v5

import blang.core.IntDistribution
import java.util.Map
import java.util.LinkedHashMap
import java.util.List
import java.util.ArrayList
import java.util.Collections
import blang.inits.experiments.Experiment
import humi.v5.DeltaMethod.Columns
import blang.inits.Arg
import humi.HumiData
import blang.types.Index
import humi.IntMixture
import blang.types.StaticUtils
import xlinear.MatrixOperations
import blang.inits.DefaultValue
import java.io.File
import blang.runtime.Runner
import briefj.BriefIO
import blang.distributions.YuleSimon
import humi.HumiStaticUtils
import blang.inits.parsing.Posix
import briefj.collections.Counter
import org.apache.commons.math3.stat.descriptive.SummaryStatistics

class M_6_Postprocess extends Experiment {
  
  @Arg
  public File execDir
  
  @Arg HumiData data
  
  @Arg   @DefaultValue("Rscript")
  public String rCmd = "Rscript"
  
  // (iteration, rna)
  val Map<Pair<Integer,Integer>,IntDistribution> samples = new LinkedHashMap
  var int nIterations = 0
  
  def load() {
    val samplesFolder = new File(execDir, Runner::SAMPLES_FOLDER)
    val rho1s = load(new File(samplesFolder, "rho1s.csv"))
    val rho2s = load(new File(samplesFolder, "rho2s.csv"))
    val pis =   load(new File(samplesFolder, "pis.csv"))
    for (sgRNA : data.targets.indices) {
      for (i : 0 ..< nIterations) {
        val key = i -> sgRNA.key
        val pi = pis.get(key)
        val rho1 = rho1s.get(key)
        val rho2 = rho2s.get(key)
        val value = IntMixture::distribution(
              StaticUtils::fixedSimplex(pi, 1.0 - pi), 
              #[ // Warning: if edited, need to change M_5_Postprocess
                YuleSimon::distribution(StaticUtils::fixedReal(rho1)), 
                YuleSimon::distribution(StaticUtils::fixedReal(rho1 + rho2))
              ]
            )
        samples.put(key, value)
      }
    }
  }
  
  def Map<Pair<Integer,Integer>,Double> load(File samples) {
    val result = new LinkedHashMap
    for (line : BriefIO::readLines(samples).indexCSV) {
      val sgRNA = Integer::parseInt(line.get("sgrna"))
      val iter = Integer::parseInt(line.get("sample"))
      nIterations = Math::max(nIterations, iter)
      val key = iter -> sgRNA
      val value = Double::parseDouble(line.get("value"))
      result.put(key, value)
    }
    return result
  }
  
  
  def mix(ArrayList<IntDistribution> distributions) {
    val unif = MatrixOperations::ones(distributions.size).div(distributions.size)
    return IntMixture::distribution(StaticUtils::fixedSimplex(unif), distributions)
  }
  
  int cutoff = 1000
  
  override run() {
    load
    
    // compute truncated means
    
    for (gene : data.genes.indices)
    for (sgRNA : data.targets.indices(gene)) {
      println(gene.key + "-" + sgRNA.key)
      // in data
      for (exp : data.experiments.indices) {
        val hist = data.histograms.get(sgRNA, exp)
        val counter = new Counter
        val counts = hist.distinctCounts
        for (int c : 1 .. cutoff) {
          if (counts.contains(c))
            counter.setCount(c, hist.frequency(c))
        }
        println("viz count empirical = " + counter.totalCount)
        counter.normalize
        println(tmean(counter))
      }
      // inferred
      val summary = new SummaryStatistics
      val cloneNumberSummary = new SummaryStatistics
      for (i : 0 ..< nIterations) {
        val dist = samples.get(Pair.of(i, sgRNA.key))
        summary.addValue(tmean(dist))
        
        //
      }
      println(summary)
    }
    
    results.flushAll
  }
  
  def double tmean(IntDistribution dist) {
    var result = 0.0
    val counter = new Counter
    for (int c : 1 .. cutoff)   
      counter.setCount(c, Math.exp(dist.logDensity(c)))
    counter.normalize
    for (int c : 1 .. cutoff) 
      result += c * counter.getCount(c)
    return result
  }
  
  def double tmean(Counter<Integer> counter) {
    var result = 0.0
    for (int c : 1 .. cutoff)
      result += c * counter.getCount(c)
    return result
  }
  
  def static void main(String [] args) {
    Experiment::start(args, Posix::parse(args), HumiStaticUtils::parsingConfigs)
  }
}