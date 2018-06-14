package corrupt.post

import corrupt.PerfectPhylo

import briefj.BriefIO
import java.io.File
import blang.inits.experiments.Experiment
import blang.inits.Arg
import blang.inits.DefaultValue

class ScoreTrees extends Experiment {
  
  @Arg File csvFile
  
  @Arg
  File referenceTree
  
  @Arg 
  @DefaultValue("value")
  String field ="value"
  
  override run() {
    parsedTreeIndicators = CLMatrixUtils::fromPhylo(PerfectPhylo::parseNewick(referenceTree)) 
    score(BriefIO.readLines(csvFile).indexCSV.map[new PerfectPhylo(it.get(field))])
  }
  
  def static void main(String [] args) {
    Experiment.startAutoExit(args)
  }
  
  var SimpleCLMatrix parsedTreeIndicators = null
  
  def void score(Iterable<PerfectPhylo> phylos) {
    var distanceOutput = results.getAutoClosedBufferedWriter("distances.csv")
    distanceOutput.append("scan,distance\n")
    var count = 0
    for (phylo : phylos) {
      
      count++
      val currentMatrix = CLMatrixUtils::fromPhylo(phylo)
      val currentDistance = CLMatrixUtils::distance(parsedTreeIndicators, currentMatrix) 
      distanceOutput.append("" + count + "," + currentDistance + "\n")
    }
  }
}