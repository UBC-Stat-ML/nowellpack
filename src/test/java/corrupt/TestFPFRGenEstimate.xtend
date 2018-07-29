package corrupt

import org.junit.Test
import org.junit.Rule
import org.junit.rules.TemporaryFolder
import blang.inits.experiments.ExperimentResults
import corrupt.post.NoiseEstimator
import java.io.File
import corrupt.post.CLMatrixUtils
import corrupt.post.NoiseEstimator.NoiseStatistics

class TestFPFRGenEstimate {
  
  @Rule
  public TemporaryFolder folder = new TemporaryFolder();
  
  @Test
  def void test() {
    val data = new GenerateData => [
      results = new ExperimentResults(folder.newFolder)
      useFPFNRates = true
      nCells = 100
      nLoci = 100
      fpRate = 0.5 
      fnRate = 0.1
    ]
    data.run
    val treeFile = new File(data.results.resultsFolder, GenerateData::TREE_FILE)
    val dataFile = new File(data.results.resultsFolder, GenerateData::DATA_FILE)
    
    val stat = new NoiseStatistics
    NoiseEstimator::estimate(PerfectPhylo::parseNewick(treeFile), CLMatrixUtils::fromCSV(dataFile), stat)
    
    println("fp=" + stat.fpRate)
    println("fn=" + stat.fnRate)
    
    
  }
}