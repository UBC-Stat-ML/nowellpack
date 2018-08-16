package corrupt.post

import corrupt.PerfectPhylo

import static extension corrupt.post.CLMatrixUtils.toCSV
import briefj.BriefIO
import java.io.File
import blang.inits.experiments.Experiment
import blang.inits.Arg
import java.util.Optional
import blang.inits.DefaultValue

class AverageCLMatrices extends Experiment {
  
  @Arg 
  public File csvFile
  
  @Arg
  public Optional<File> referenceTree = Optional.empty
  
  @Arg                @DefaultValue("true")
  public boolean logisticTransform = true;
  
  @Arg 
  @DefaultValue("value")
  public String field ="value"
  
  public static val OUTPUT_NAME = "average.csv"
  override run() {
    if (referenceTree.present)
      parsedTreeIndicators = CLMatrixUtils::fromPhylo(PerfectPhylo::parseNewick(referenceTree.get)) 
    averageTipIndicators(BriefIO.readLines(csvFile).indexCSV.map[new PerfectPhylo(it.get(field))], logisticTransform)
    result.toCSV(results.getFileInResultFolder(OUTPUT_NAME), parsedTreeIndicators) 
  }
  
  def static void main(String [] args) {
    Experiment.startAutoExit(args)
  }
  
  var SimpleCLMatrix parsedTreeIndicators = null
  var SimpleCLMatrix result = null
  
  /**
   * Null if empty.
   */
  def void averageTipIndicators(Iterable<PerfectPhylo> phylos, boolean logisticTransform) {
    var distanceOutput = 
      if (referenceTree === null) 
        null 
      else 
        results.getAutoClosedBufferedWriter("distances.csv")
    if (parsedTreeIndicators !== null)
      distanceOutput.append("scan,distance\n")
    var count = 0
    var finalDistance = Double.NaN
    for (phylo : phylos) {
      count++
      if (result === null) {
        result = new SimpleCLMatrix(phylo.cells, phylo.loci)
        if (parsedTreeIndicators !== null)
          distanceOutput.append("" + count + "," + distance(parsedTreeIndicators, result, count) + "\n")
      }
      result += phylo
      if (parsedTreeIndicators !== null) {
        finalDistance = distance(parsedTreeIndicators, result, count)
        distanceOutput.append("" + count + "," + finalDistance + "\n")
      }
    }
    if (!Double.isNaN(finalDistance)) {
      println("distance = " + finalDistance)
    }
    if (result !== null) {
      result /= count
      if (logisticTransform)
        result.logisticTransform
    }
  }
  
  private def double distance(SimpleCLMatrix refTree, SimpleCLMatrix sum, double count) {
    CLMatrixUtils::checkCompatible(refTree, sum)
    return CLMatrixUtils::distance(refTree.matrix, sum.matrix/count) 
  }
}