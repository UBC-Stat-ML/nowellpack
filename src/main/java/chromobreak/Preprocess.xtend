package chromobreak

import blang.inits.experiments.Experiment
import blang.inits.Arg
import java.io.File
import blang.inits.DefaultValue

class Preprocess extends Experiment {
  
  @Arg File reads
  @Arg File gc
  
  @Arg            @DefaultValue("INF")
  public int maxNCells = Integer.MAX_VALUE
  
  override run() {
    processReads
    processGC
  }
  
  def processReads() {
    val tidifyFolder = results.child("tidyReads")
    val configs = experimentConfigs
    val tidify = new TidifyCounts => [
      cellSpecificFiles = true
      countFile = reads
      useInteger = false
      results = tidifyFolder
      experimentConfigs = configs
      nCells = maxNCells
    ]
    tidify.run
  }
  
  def processGC() {
    val tidifyFolder = results.child("tidyGC")
    val configs = experimentConfigs
    val tidify = new TidifyCounts => [
      cellSpecificFiles = false
      countFile = gc
      useInteger = false
      results = tidifyFolder
      experimentConfigs = configs
      nCells = 1
    ]
    tidify.run
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}