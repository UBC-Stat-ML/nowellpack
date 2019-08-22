package humi.freq

import bayonet.smc.ParticlePopulation
import bayonet.smc.ResamplingScheme
import binc.Command
import blang.inits.Arg
import blang.inits.DefaultValue
import blang.inits.experiments.Experiment
import blang.inits.parsing.Posix
import blang.types.Plated
import briefj.BriefIO
import humi.CountFrequencies
import humi.HumiData
import humi.HumiStaticUtils
import humi.SimpleCountFrequencies
import java.util.ArrayList
import java.util.List
import java.util.Random
import humi.Preprocess

class GenerateData extends Experiment {
  @Arg HumiData data
  @Arg List<String> experiments
  @Arg List<Double> lambdas
  @Arg Plated<Integer> initialPopCounts  // TODO: resample this as well!
  
  @Arg        @DefaultValue("1")
  Random random = new Random(1)
  
  @Arg              @DefaultValue("false")
  boolean justComputeLambdaBound = false
  
  @Arg   @DefaultValue("Rscript")
  public String rCmd = "Rscript"
  
  @Arg              @DefaultValue("1.0")
  public double relativeDataSize = 1.0
  
  override run() {
    if (experiments.size !== lambdas.size)
      throw new RuntimeException
    for (experiment : data.experiments.indices) {
      val lambda = lambda(experiment.key)
      if (lambda !== null) { 
        val expIndex = experiment.plate.name -> experiment.key
        var bound = 0.0
        for (gene : data.genes.indices)
        for (sgRNA : data.targets.indices(gene)) {
          val rnaIndex = sgRNA.plate.name -> sgRNA.key
          val geneIndex = gene.plate.name -> gene.key
          val trueHistogram = data.histograms.get(sgRNA, experiment)
          val nNonZeroClones = trueHistogram.nDataPoints
          val curBound = nNonZeroClones as double / initialPopCounts.get(sgRNA)
          if (justComputeLambdaBound)
            println('''bound = S_t / J_t = «trueHistogram.nDataPoints» / «initialPopCounts.get(sgRNA)» = «curBound»''')
          bound = Math::max(bound, curBound )
          if (!justComputeLambdaBound) {
            val double imputedZeros = lambda * initialPopCounts.get(sgRNA) - trueHistogram.nDataPoints
            if (justComputeLambdaBound)
              println("\tImputed zeros " + imputedZeros)
            if (imputedZeros < 0) 
              throw new RuntimeException('''Try lower lambda. At sgRNA «sgRNA.key», experiment «experiment.key» expected non-negative value: «lambda * initialPopCounts.get(sgRNA)» - «trueHistogram.nDataPoints» = «imputedZeros»''')
            val resampledHistogram = resample(trueHistogram, imputedZeros)
            record(sgRNA.key, trueHistogram, imputedZeros, resampledHistogram)
            results.getTabularWriter("complete").write(expIndex, rnaIndex, geneIndex, Preprocess::HISTOGRAM_COLUMN -> resampledHistogram)
            resampledHistogram.dropZeros
            if (justComputeLambdaBound)
              println("\tnNonZeroClones before and after resampling: " + nNonZeroClones + " " + resampledHistogram.nDataPoints)
            results.getTabularWriter("censored").write(expIndex, rnaIndex, geneIndex, Preprocess::HISTOGRAM_COLUMN -> resampledHistogram)
          }
        }
        if (justComputeLambdaBound)
          println("Experiment " + experiment.key + " bound: " + bound)
        }
    }
    if (!justComputeLambdaBound)
      plot
  }
  
  def Double lambda(String string) {
    for (i : 0 ..< experiments.size) {
      if (experiments.get(i) == string) return lambdas.get(i)
    }
    return null
  }
  
  def plot() {
    results.flushAll
    val plotResults = results.child("plots")
    val scriptFile = plotResults.getFileInResultFolder("script.r")
    BriefIO::write(scriptFile, '''
      require("ggplot2")
      data <- read.csv("«results.getFileInResultFolder("all.csv").absolutePath»")
      p <- ggplot(data, aes(x = cloneSize, y = log10(freq), shape = factor(type), colour = factor(type))) + 
        geom_point() + 
        facet_grid(~ sgRNA)
      ggsave("«plotResults.getFileInResultFolder("freqs.pdf").absolutePath»", width = 3000, limitsize = F)
    ''')
    Command.call(Command.cmd(rCmd).appendArg(scriptFile.getAbsolutePath()))
  }
  
  def record(Integer id, CountFrequencies trueHistogram, double imputedZeros, SimpleCountFrequencies resampledHistogram) {
    // record true
    val writer = results.getTabularWriter("all").child("sgRNA", id)
    val originals = writer.child("type", "original")
    for (count : trueHistogram.distinctCounts)
      originals.write(
        "cloneSize" -> count,
        "freq" -> trueHistogram.frequency(count)
      )
    originals.write(
      "cloneSize" -> 0,
      "freq" -> imputedZeros
    )
    for (count : resampledHistogram.distinctCounts)
      writer.write(
        "type" -> "resampled",
        "cloneSize" -> count,
        "freq" -> resampledHistogram.frequency(count)
      )
    results.flushAll
  }
  
  def resample(CountFrequencies frequencies, double imputedZeros) {
    val logWeights = newDoubleArrayOfSize(frequencies.distinctCounts.size + 1)
    val particles = new ArrayList<Integer>
    logWeights.set(0, Math::log(imputedZeros))
    particles.add(0)
    var i = 1
    for (count : frequencies.distinctCounts) {
      logWeights.set(i++, Math::log(frequencies.frequency(count)))
      particles.add(count)
    }
    var nTimesToResample = (imputedZeros + frequencies.nDataPoints) * relativeDataSize
    val preResampling = ParticlePopulation.buildDestructivelyFromLogWeights(logWeights, particles, 0.0)
    val resampled = ResamplingScheme.MULTINOMIAL.resample(random, preResampling.normalizedWeights, particles, nTimesToResample.intValue);
    return new SimpleCountFrequencies => [addAll(resampled)]
  }
  
  def double [] getNormalizedWeights(ParticlePopulation<?> pp) {
    val result = newDoubleArrayOfSize(pp.nParticles)
    for (i : 0 ..< pp.nParticles)
      result.set(i, pp.getNormalizedWeight(i))
    return result
  }
  
  def static void main(String [] args) {
    Experiment::start(args, Posix::parse(args), HumiStaticUtils::parsingConfigs)
  }
}