package corrupt.distances

import java.io.File
import blang.inits.Arg
import blang.inits.experiments.Experiment
import corrupt.PerfectPhylo
import blang.inits.DefaultValue
import briefj.BriefIO
import briefj.collections.Counter

class L1Decoder extends Experiment {
  @Arg
  public File samples
  
  @Arg 
  @DefaultValue("value")
  String field ="value"
  
  @Arg  @DefaultValue("LocusEdgeStat")
  public TreeStatistic<?> stat = new LocusEdgeStat
  
  def static <T> PerfectPhylo run(
    Iterable<PerfectPhylo> candidates, 
    TreeStatistic<T> statistic,
    Counter<T> posteriorStatistics
  ) {
    for (posterior : posteriorStatistics.entries.values)
      if (posterior < 0.0 || posterior > 1.0)
        throw new RuntimeException
    var PerfectPhylo result = null
    var minLoss = Double.POSITIVE_INFINITY
    for (candidate : candidates) {
      val stats = statistic.compute(candidate)
      var sum = 0.0
      for (stat : stats) {
        val posterior = posteriorStatistics.getCount(stat)
        sum += 1.0 - 2.0 * posterior
      }
      if (sum < minLoss) {
        minLoss = sum
        result = candidate
      }
    }
    return result
  }
  
  override run() {
    val posterior = posterior(trees, stat)
    val result = run(trees, stat as TreeStatistic, posterior as Counter)
    BriefIO::write(results.getFileInResultFolder("consensus.newick"), result.toNewick) 
  }
  
  def static <T> Counter<T> posterior(Iterable<PerfectPhylo> trees, TreeStatistic<T> stat) {
    val result = new Counter
    var nTrees = 0
    for (tree : trees) {
      result.incrementAll(stat.compute(tree), 1.0)
      nTrees++
    }
    for (key : result.entries.keySet) {
      val cur = result.getCount(key)
      result.setCount(key, cur / nTrees)
    }
    return result
  }
  
  def trees() {
    return BriefIO::readLines(samples).indexCSV.map[new PerfectPhylo(it.get(field))]
  }
  
  def public static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}