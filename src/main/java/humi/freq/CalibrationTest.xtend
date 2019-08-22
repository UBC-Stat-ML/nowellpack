package humi.freq

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
  
  @Arg File intervals
  
  @Arg @DefaultValue("sgrna")
  String sgRNAField = "sgrna"
  
  override run() {
    val truth = loadReference()
    val widths = new SummaryStatistics
    val coverage = new SummaryStatistics
    for (line : BriefIO::readLines(intervals).indexCSV) { 
      val sgRNA = Integer::parseInt(line.get(sgRNAField))
      if (!truth.containsKey(sgRNA)) throw new RuntimeException
      val left = Double::parseDouble(line.get(DeltaMethod.Columns.logRatioLeftBound.toString))
      val right = Double::parseDouble(line.get(DeltaMethod.Columns.logRatioRightBound.toString))
      widths.addValue(right - left)
      val curRef = truth.get(sgRNA)
      coverage.addValue(if (curRef >= left && curRef <= right) 1.0 else 0.0)
      truth.remove(sgRNA) 
    }
    println("Average width: " + widths.mean)
    println("Actual coverage: " + coverage.mean)
  }
  
  def loadReference() {
    val result = new LinkedHashMap<Integer, Double>
    for (line : BriefIO::readLines(reference).indexCSV) { 
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