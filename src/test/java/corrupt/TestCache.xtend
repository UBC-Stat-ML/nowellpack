package corrupt

import org.junit.Test
import blang.runtime.Runner
import bayonet.distributions.Random
import corrupt.post.ReadOnlyCLMatrix
import corrupt.post.CLMatrixUtils
import blang.engines.internals.factories.SCM
import blang.inits.experiments.ExperimentResults
import org.junit.rules.TemporaryFolder
import org.junit.Rule
import corrupt.post.NoisyBinaryCLMatrix

class TestCache {
  @Rule
  public TemporaryFolder folder = new TemporaryFolder();
  
  @Test def void testReadOnlyMtxCache() {
    val nCells = 10
    val nLoci = 20
    val rand = new Random(1)
    CorruptPhylo::testCacheCorrectness = true
    // generate data
    val phylo = PerfectPhylo::generateUniform(nCells, nLoci, rand)
    val matrix = ReadOnlyCLMatrix::readOnly(CLMatrixUtils::syntheticInclusionPrs(rand, phylo, 0.5))   
    // create inference machinery
    val modelBuilder = new FixedMatrixModel.Builder().setTipInclusionProbabilities(matrix)
    val runner = new Runner(modelBuilder)
    runner.results = new ExperimentResults(folder.newFolder) 
    runner.engine = new SCM => [
      nParticles = 20
      results = runner.results
    ]
    runner.run
    CorruptPhylo::testCacheCorrectness = false 
  }
  
  @Test def void testNoisyMtxCache() {
    val nCells = 10
    val nLoci = 20
    val rand = new Random(1)
    CorruptPhylo::testCacheCorrectness = true
    // generate data
    val phylo = PerfectPhylo::generateUniform(nCells, nLoci, rand)
    val binaryMatrix = ReadOnlyCLMatrix::readOnly(CLMatrixUtils::syntheticPerturbedBinaryMatrix(rand, phylo, 0.2, 0.3)) 
    val noisy = new NoisyBinaryCLMatrix(binaryMatrix, 0.1, 0.05) 
    // create inference machinery
    val modelBuilder = new FixedMatrixModel.Builder().setTipInclusionProbabilities(noisy)
    val runner = new Runner(modelBuilder)
    runner.results = new ExperimentResults(folder.newFolder) 
    runner.engine = new SCM => [
      nParticles = 20
      results = runner.results
    ]
    runner.run
    CorruptPhylo::testCacheCorrectness = false 
  }
}
