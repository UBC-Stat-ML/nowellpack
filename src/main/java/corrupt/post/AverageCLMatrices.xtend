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
  
  @Arg File csvFile
  
  @Arg
  Optional<File> referenceTree
  
  @Arg 
  @DefaultValue("value")
  String field ="value"
  
  public static val OUTPUT_NAME = "average.csv"
  override run() {
    if (referenceTree.present)
      parsedTreeIndicators = CLMatrixUtils::fromPhylo(PerfectPhylo::parseNewick(referenceTree.get)) 
    averageTipIndicators(BriefIO.readLines(csvFile).indexCSV.map[new PerfectPhylo(it.get(field))])
    result.toCSV(results.getFileInResultFolder(OUTPUT_NAME), parsedTreeIndicators) 
  }
  
  def static void main(String [] args) {
    Experiment.start(args)
  }
  
  var SimpleCLMatrix parsedTreeIndicators = null
  var SimpleCLMatrix result = null
  
  /**
   * Null if empty.
   */
  def void averageTipIndicators(Iterable<PerfectPhylo> phylos) {
    var distanceOutput = 
      if (referenceTree === null) 
        null 
      else 
        results.getAutoClosedBufferedWriter("distances.csv") => [ println("distance") ]
    var count = 0
    var finalDistance = Double.NaN
    for (phylo : phylos) {
      count++
      if (result === null) {
        result = new SimpleCLMatrix(phylo.cells, phylo.loci)
        if (parsedTreeIndicators !== null)
          distanceOutput.append(distance(parsedTreeIndicators, result, count) + "\n")
      }
      result += phylo
      if (parsedTreeIndicators !== null) {
        finalDistance = distance(parsedTreeIndicators, result, count)
        distanceOutput.append(finalDistance + "\n")
      }
    }
    if (!Double.isNaN(finalDistance)) {
      println("distance = " + finalDistance)
    }
    if (result !== null)
      result /= count
  }
  
  private def double distance(SimpleCLMatrix refTree, SimpleCLMatrix sum, double count) {
    CLMatrixUtils::checkCompatible(refTree, sum)
    return CLMatrixUtils::distance(refTree.matrix, sum.matrix/count) 
  }
  
  
}