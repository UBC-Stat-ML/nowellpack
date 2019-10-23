package humi

import blang.runtime.internals.DefaultPostProcessor
import java.io.File
import blang.runtime.Runner
import java.util.List
import java.util.Map
import java.util.LinkedHashMap
import briefj.BriefIO
import blang.inits.experiments.tabwriters.TidySerializer
import briefj.BriefMaps
import blang.inits.Arg
import java.util.ArrayList
import org.apache.commons.math3.stat.descriptive.SummaryStatistics
import blang.types.Index
import humi.freq.DeltaMethod.Columns
import java.util.Collections
import blang.inits.DefaultValue
import blang.inits.experiments.Experiment
import blang.inits.parsing.Posix

class HumiPostProcessor extends DefaultPostProcessor {
  
  @Arg public HumiData data
  
  @Arg                @DefaultValue("0.9")
  public double credibleIntervalPr = 0.9
  
  var nIterations = 0
  
  override run() {
    super.run
    computeIntervals
    results.flushAll
    // GoF diagnostic summary
    for (stat : GofStat.values)
      gofSummary(stat)
    // intervals
    HumiStaticUtils::plotIntervals(results, data, rCmd, false, "Bayesian hierarchical model credible intervals")
    // TODO: poset
  }
  
  // Note: currently assuming by default it is non-experiment-specific, except for hard coded exception for visibleCloneNumbers (*)
  static enum GofStat { truncatedMeans, truncatedSqMeans, visibleCloneNumbers }
  
  def void gofSummary(GofStat stat) {
    val samplesDir = new File(blangExecutionDirectory.get, Runner::SAMPLES_FOLDER)
    
    val sampleFile = new File(samplesDir, stat.toString + ".csv")
    val samples = new LinkedHashMap<Object,List<Double>>  // sgrna, experiment -> samples
    for (line : BriefIO::readLines(sampleFile).indexCSV) {
      val sgRNAName = data.targets.name.toString 
      val sgRNAIdString = line.get(sgRNAName)
      val sgRNA = Integer::parseInt(sgRNAIdString) 
      val expString = line.get(data.experiments.name.toString)
      val key = if (stat == GofStat.visibleCloneNumbers) Pair.of(sgRNA, expString) else sgRNA  
      val value = Double::parseDouble(line.get(TidySerializer::VALUE))
      BriefMaps::getOrPutList(samples, key).add(value)
    }
    
    for (Index<String> experiment : data.experiments.indices) {
      val coverageSummary = new SummaryStatistics
      val widthSummary = new SummaryStatistics
      for (Index<Integer> target : data.targets.indices) {
        val key = if (stat == GofStat.visibleCloneNumbers) Pair.of(target.key, experiment.key) else target.key
        val values = samples.get(key)
        if (values !== null) { // can be null when inference was performed on a subset and looking at visibleCloneNumbers
          Collections.sort(values)
          val a = 1.0 - credibleIntervalPr
          val leftBound = values.get(((a/2.0) * values.size) as int)
          val rightBound = values.get(((1.0 - (a/2)) * values.size) as int)
          
          val observedHist = data.histograms.get(target, experiment)
          val observed = 
            if (stat == GofStat.truncatedMeans) DistributionSummary::mean(DistributionSummary::truncatedNormalizedCounter(observedHist))
            else if (stat == GofStat.truncatedSqMeans) DistributionSummary::meanSq(DistributionSummary::truncatedNormalizedCounter(observedHist))
            else if (stat == GofStat.visibleCloneNumbers) observedHist.nDataPoints
            else throw new RuntimeException // see (*) above too
          val contains = leftBound <= observed && observed <= rightBound
          coverageSummary.addValue(if (contains) 1.0 else 0.0)
          widthSummary.addValue(rightBound - leftBound)
        }
      }
      if (coverageSummary.n > 0)
        results.getTabularWriter("gof").write(
          "referenceDataset" -> experiment.key,
          "gofStatistic" -> stat.toString,
          "width" -> widthSummary.mean,
          "theoreticalCoverage" -> credibleIntervalPr,
          "actualCoverage" -> coverageSummary.mean
        )
    }
  }
  
  def void computeIntervals() {
    val samplesDir = new File(blangExecutionDirectory.get, Runner::SAMPLES_FOLDER)
    val condWinMeansFile = new File(samplesDir, "conditionalWinsorizedMeans.csv")
    val samples = read(condWinMeansFile)
    val controlWMeans = controlWMeans(samples)
    for (gene : data.genes.indices)
      for (sgRNA : data.targets.indices(gene)) 
        if (samples.containsKey(sgRNA.key)) {
          val List<Double> wMeans = samples.get(sgRNA.key)
          val List<Double> logRatios = new ArrayList
          for (i : 0 ..< controlWMeans.size)
            logRatios.add(Math.log(wMeans.get(i) / controlWMeans.get(i)))
          credibleIntervals(sgRNA, gene, logRatios)
        }
  }
  
  def credibleIntervals(Index<Integer> sgRNA, Index<String> gene, List<Double> values) {
    Collections.sort(values)
    if (credibleIntervalPr < 0.0 || credibleIntervalPr > 1.0) 
      throw new RuntimeException
    val a = 1.0 - credibleIntervalPr
    val int left = ((a/2.0) * values.size) as int
    val int middle = (0.5 * values.size) as int
    val int right = ((1.0 - (a/2)) * values.size) as int
    results.getTabularWriter("estimates").write(
      sgRNA.plate.name -> sgRNA.key,
      gene.plate.name -> gene.key,
      Columns::logRatioLeftBound -> values.get(left),
      Columns::logRatioRightBound -> values.get(right),
      Columns::logRatio -> values.get(middle)
    )
  }
  
  def List<Double> controlWMeans(Map<Integer,List<Double>> samples) {
    val result = new ArrayList<Double>
    for (iter : 0 ..< nIterations) {
      val stats = new SummaryStatistics
      for (gene : data.genes.indices) {
        if (data.isControl(gene)) {
          for (sgRNA : data.targets.indices(gene)) 
            if (samples.containsKey(sgRNA.key)) { //  need to check in case inference was on subset
              val list = samples.get(sgRNA.key)
              if (list.size !== nIterations)
                throw new RuntimeException
              stats.addValue(list.get(iter))
            }
        }
      }
      result.add(stats.mean)
    }
    return result
  }
  
  def Map<Integer,List<Double>> read(File csvFile) {
    val result = new LinkedHashMap
    for (line : BriefIO::readLines(csvFile).indexCSV) {
      val sgRNAName = data.targets.name.toString 
      val sgRNAIdString = line.get(sgRNAName)
      val sgRNA = Integer::parseInt(sgRNAIdString) 
      val value = Double::parseDouble(line.get(TidySerializer::VALUE))
      val list = BriefMaps::getOrPutList(result, sgRNA)
      list.add(value)
      nIterations = Math.max(nIterations, list.size)
    }
    return result
  }
  
  def static void main(String [] args) {
    Experiment::start(args, Posix::parse(args), HumiStaticUtils::parsingConfigs)
  }
}