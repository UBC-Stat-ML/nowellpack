package blang.inits.experiments

import org.junit.Test
import bayonet.distributions.Random
import corrupt.post.ReadOnlyCLMatrix
import corrupt.post.CLMatrixUtils
import corrupt.post.NoisyBinaryCLMatrix
import blang.runtime.Runner
import org.junit.Rule
import org.junit.rules.TemporaryFolder
import blang.mcmc.internals.SamplerBuilderOptions
import java.io.File
import org.apache.commons.math3.stat.descriptive.SummaryStatistics
import briefj.BriefIO
import blang.types.StaticUtils
import blang.engines.internals.factories.PT
import java.util.Optional
import corrupt.CorruptGibbsSampler
import corrupt.GlobalErrorModel
import corrupt.PerfectPhylo
import corrupt.CorruptPhylo
import org.junit.Assert

class ParamEstimation {
  @Rule
  public TemporaryFolder folder = new TemporaryFolder();
  
  @Test
  def void test() {
    val truth = #{"fpr" -> 0.15, "fnr" -> 0.17}
    val nCells = 1000
    val nLoci = 100 
    val rand = new Random(1)
    // generate data
    val phylo = PerfectPhylo::generateUniform(nCells, nLoci, rand)
    val binaryMatrix = ReadOnlyCLMatrix::readOnly(CLMatrixUtils::syntheticPerturbedBinaryMatrix(rand, phylo, truth.get("fpr"), truth.get("fnr"))) 
    val fpr = StaticUtils::latentReal
    fpr.set(0.1)
    val fnr = StaticUtils::latentReal
    fnr.set(0.1)
    val noisy = new NoisyBinaryCLMatrix(binaryMatrix, fpr, fnr)
    // create inference machinery
    val modelBuilder = new GlobalErrorModel.Builder().setBinaryMatrix(binaryMatrix).setFpr(fpr).setFnr(fnr).setPhylo(new CorruptPhylo(phylo, noisy))
    val runner = new Runner(modelBuilder)
    runner.stripped = true
    runner.samplers = new SamplerBuilderOptions => [
      excluded.add(CorruptGibbsSampler) 
    ]
    runner.results = new ExperimentResults(folder.newFolder) 
    runner.engine = new PT => [
      results = runner.results 
      nScans = 1000
      nChains = Optional.of(1)
    ]
    runner.run
    runner.results.closeAll
    for (stat : truth.keySet) {
      val samplesFile = new File(runner.results.resultsFolder, "samples/" + stat + ".csv")
      val summStat = new SummaryStatistics
      for (line : BriefIO.readLines(samplesFile).indexCSV)
        summStat.addValue(Double.parseDouble(line.get("value")))
      println(stat + " = " + summStat.mean)
      Assert.assertEquals(truth.get(stat), summStat.mean, 0.01) 
    }
  }
}