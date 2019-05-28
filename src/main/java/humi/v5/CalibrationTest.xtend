package humi.v5

import blang.inits.experiments.Experiment
import blang.inits.parsing.Posix
import humi.HumiStaticUtils
import blang.inits.Arg
import blang.inits.DefaultValue
import java.io.File
import briefj.BriefIO
import java.util.LinkedHashMap
import org.apache.commons.math3.stat.descriptive.SummaryStatistics

class CalibrationTest  extends Experiment {
  
  @Arg File reference
  @Arg String referenceDataset
  
  @Arg File intervals
  
  @Arg @DefaultValue("dataset")
  String datasetField = "dataset"
  
  @Arg @DefaultValue("sgrna")
  String sgRNAField = "sgrna"
  
  override run() {
    val truth = loadReference()
    val radia = new SummaryStatistics
    val coverage = new SummaryStatistics
    for (line : BriefIO::readLines(intervals).indexCSV) { //}.filter[get(datasetField) != referenceDataset]) {
      val sgRNA = Integer::parseInt(line.get(sgRNAField))
      val logRatio = Double::parseDouble(line.get(DeltaMethod.Columns.logRatio.toString))
      val radius = Double::parseDouble(line.get(DeltaMethod.Columns.logRatioIntervalRadius.toString))
      radia.addValue(radius)
      coverage.addValue(if (Math::abs(logRatio - truth.get(sgRNA)) < radius) 1.0 else 0.0)
    }
    println(radia.mean)
    println(coverage.mean)
  }
  
  def loadReference() {
    val result = new LinkedHashMap<Integer, Double>
    for (line : BriefIO::readLines(reference).indexCSV) { //}.filter[get(datasetField) == referenceDataset]) {
      val sgRNA = Integer::parseInt(line.get(sgRNAField))
      val logRatio = Double::parseDouble(line.get(DeltaMethod.Columns.logRatio.toString))
      result.put(sgRNA, logRatio)
    }
    return result
  }
  
  def static void main(String [] args) {
    Experiment::start(args, Posix::parse(args), HumiStaticUtils::parsingConfigs)
  }
  
}