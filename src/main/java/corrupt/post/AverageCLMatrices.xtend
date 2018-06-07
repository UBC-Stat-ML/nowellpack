package corrupt.post

import corrupt.PerfectPhylo
import briefj.BriefIO
import java.io.File
import blang.inits.experiments.Experiment
import blang.inits.Arg
import java.util.Optional

class AverageCLMatrices extends Experiment {
  
  @Arg File csvFile
  @Arg Optional<String> field
  
  public static val OUTPUT_NAME = "average.csv"
  override run() {
    averageTipIndicators(csvFile, field.orElse(null)).toCSV(results.getFileInResultFolder(OUTPUT_NAME))
  }
  
  def static void main(String [] args) {
    Experiment.start(args)
  }
  
  /**
   * Null if empty.
   */
  def static SimpleCLMatrix averageTipIndicators(Iterable<PerfectPhylo> phylos) {
    var SimpleCLMatrix result = null
    var count = 0
    for (phylo : phylos) {
      count++
      if (result === null) 
        result = new SimpleCLMatrix(phylo.cells, phylo.loci)
      result += phylo
    }
    if (result === null)
      return null
    result /= count
    return result
  }
  
  def static SimpleCLMatrix averageTipIndicators(File csvFile) {
    return averageTipIndicators(csvFile, null) 
  }
  def static SimpleCLMatrix averageTipIndicators(File csvFile, String _field) {
    val field = if (_field === null) "value" else _field
    return averageTipIndicators(BriefIO.readLines(csvFile).indexCSV.map[new PerfectPhylo(it.get(field))])
  }
}