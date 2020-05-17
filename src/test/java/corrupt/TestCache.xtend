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
import java.util.LinkedHashSet
import java.util.Set

class TestCache {
  @Rule
  public TemporaryFolder folder = new TemporaryFolder();
  
  @Test def void testConcatenationCache() {
    val globalParameterization = true
    val nCells = 10
    val nLoci = 20
    val rand = new Random(1)
    CorruptPhylo::testCacheCorrectness = true
    // generate tree
    val phylo = PerfectPhylo::generateUniform(nCells, nLoci, rand)
    // split loci
    val Set<Locus> set1 = new LinkedHashSet
    val Set<Locus> set2 = new LinkedHashSet
    var i = 0
    for (locus : phylo.loci)
      (if (i++ < 5) set1 else set2).add(locus)
    // generate data (read only)
    
    val snvs = ReadOnlyCLMatrix::readOnly(CLMatrixUtils::syntheticInclusionPrs(rand, phylo, 0.5, set1)) 
    // generate data (noisy)
    val binaryMatrix = CLMatrixUtils::syntheticPerturbedBinaryMatrix(rand, phylo, 0.2, 0.3, set2)
    val modelBuilder = new ConcatenationModel.Builder()
      .setCopyNumberMarkers(binaryMatrix)
      .setSingleNucleotideMarkers(snvs)
      .setGlobalParameterization(globalParameterization)
      .setFprBound(0.2)
      .setFnrBound(0.2)
      .setFpr(corrupt.CorruptStaticUtils.initializedLatentErrors(0.0, 0.2, nLoci, globalParameterization))
      .setFnr(corrupt.CorruptStaticUtils.initializedLatentErrors(0.0, 0.2, nLoci, globalParameterization))
    val runner = new Runner(modelBuilder)
    runner.results = new ExperimentResults(folder.newFolder) 
    runner.engine = new SCM => [
      nParticles = 20
      results = runner.results
    ]
    runner.run
    CorruptPhylo::testCacheCorrectness = false 
  }
  
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
    val binaryMatrix = CLMatrixUtils::syntheticPerturbedBinaryMatrix(rand, phylo, 0.2, 0.3)
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
