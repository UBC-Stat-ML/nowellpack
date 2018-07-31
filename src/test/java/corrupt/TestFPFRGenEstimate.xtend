package corrupt

import org.junit.Test
import org.junit.Rule
import org.junit.rules.TemporaryFolder
import blang.inits.experiments.ExperimentResults
import corrupt.post.NoiseEstimator
import java.io.File
import corrupt.post.CLMatrixUtils
import corrupt.post.NoiseStatistics
import org.junit.Assert

class TestFPFRGenEstimate {
  
  @Rule
  public TemporaryFolder folder = new TemporaryFolder();
  
  @Test
  def void test() {
    for (fp : #[0.0, 0.2, 0.5]) {
      for (fn : #[0.0, 0.17, 0.9]) {
        val data = new GenerateData => [
          results = new ExperimentResults(folder.newFolder)
          useFPFNRates = true
          nCells = 100
          nLoci = 100
          fpRate = fp 
          fnRate = fn
        ]
        data.run
        val treeFile = new File(data.results.resultsFolder, GenerateData::TREE_FILE)
        val dataFile = new File(data.results.resultsFolder, GenerateData::DATA_FILE)
        
        val stat = new NoiseStatistics
        NoiseEstimator::estimate(PerfectPhylo::parseNewick(treeFile), CLMatrixUtils::fromCSV(dataFile), stat)
        Assert.assertEquals(fp, stat.fpRate, 0.01)
        Assert.assertEquals(fn, stat.fnRate, 0.01)
      }
    }
  }
}