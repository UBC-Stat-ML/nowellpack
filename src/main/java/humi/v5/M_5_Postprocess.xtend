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

class M_5_Postprocess extends Experiment {
  
  @Arg                  @DefaultValue("0.95")
  public double winsorizedTailCutoff = 0.95
  
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
  
  def computeIntervals() {
    val List<Double> controlWMeans = controlWMeans()
    for (gene : data.genes.indices)
      for (sgRNA : data.targets.indices(gene)) {
        val List<Double> wMeans = wMeans(sgRNA.key)
        val List<Double> logRatios = new ArrayList
        for (i : 0 ..< controlWMeans.size)
          logRatios.add(Math.log(wMeans.get(i) / controlWMeans.get(i)))
        credibleIntervals(sgRNA, gene, logRatios)
      }
  }
  
  def credibleIntervals(Index<Integer> sgRNA, Index<String> gene, List<Double> values) {
    Collections.sort(values)
    val int left = (0.025 * values.size) as int
    val int middle = (0.5 * values.size) as int
    val int right = (0.975 * values.size) as int
    results.getTabularWriter("estimates").write(
      sgRNA.plate.name -> sgRNA.key,
      gene.plate.name -> gene.key,
      Columns::logRatioLeftBound -> values.get(left),
      Columns::logRatioRightBound -> values.get(right),
      Columns::logRatio -> values.get(middle)
    )
  }
    
  def controlWMeans() {
    val result = new ArrayList<Double>
    for (i : 0 ..< nIterations) {
      // build the mixture
      val dists = new ArrayList<IntDistribution>
      for (gene : data.genes.indices) {
        if (data.isControl(gene)) {
          for (sgRNA : data.targets.indices(gene)) {
            val key = i -> sgRNA.key
            val dist = samples.get(key)
            dists.add(dist)
          }
        }
      }
      result.add(conditionalWinsorizedMean(mix(dists), winsorizedTailCutoff))
    }
    return result
  }
  
  def mix(ArrayList<IntDistribution> distributions) {
    val unif = MatrixOperations::ones(distributions.size).div(distributions.size)
    return IntMixture::distribution(StaticUtils::fixedSimplex(unif), distributions)
  }
  
  def List<Double> wMeans(Integer sgRNA) {
    val result = new ArrayList
    for (i : 0 ..< nIterations) {
      val key = i -> sgRNA
      val dist = samples.get(key)
      result.add(conditionalWinsorizedMean(dist, winsorizedTailCutoff))
    }
    return result
  }
  
  def double conditionalWinsorizedMean(IntDistribution distribution, double p) {
    if (p < 0.5 || p > 1.0) throw new RuntimeException
    var sum = 0.0
    var x = 0
    val normalization = 1.0 - Math::exp(distribution.logDensity(0))
    while (sum/normalization < p) {
      x++  
      sum += Math::exp(distribution.logDensity(x))
    }
    val cutOff = x
    var result = 0.0
    var mass = 0.0
    for (y : 0 ..< cutOff) {
      val currentMass = Math::exp(distribution.logDensity(y))
      mass += currentMass
      result += y * currentMass
    }
    result += cutOff * (1.0 - mass)
    return result
  }
  
  override run() {
    load
    computeIntervals
    results.flushAll
    DeltaMethod::plot(results, data, rCmd, false, "Bayesian hierarchical model credible intervals")
  }
  
  def static void main(String [] args) {
    Experiment::start(args, Posix::parse(args), HumiStaticUtils::parsingConfigs)
  }
}