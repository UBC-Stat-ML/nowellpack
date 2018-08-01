package corrupt.post

import blang.inits.experiments.Experiment
import blang.inits.Arg
import corrupt.PerfectPhylo
import java.util.ArrayList
import java.util.Collections
import corrupt.Locus
import java.util.Comparator
import corrupt.Greedy
import blang.inits.DefaultValue

class RankLoci extends Experiment {
  @Arg PerfectPhylo phylo
  @Arg ReadOnlyCLMatrix binaryMatrix
  
  @Arg         @DefaultValue("20")
  public int thinningPeriod = 20
  
  def stat(Locus locus) {
    val stat = new NoiseStatistics
    stat.add(locus, phylo.getTips(locus), binaryMatrix)
    return stat
  }
  
  override run() {
    val sortedLoci = new ArrayList(phylo.loci)
    Collections::sort(
      sortedLoci, 
      Comparator::comparing[Locus locus | stat(locus).fpRate]
    )
    val incrementalMatrix = new SimpleCLMatrix(phylo.cells, phylo.loci)
    var iteration = 0
    val treesWriter = results.getTabularWriter("trees")
    for (locus : sortedLoci) {
      val tips = phylo.getTips(locus)
      for (entry : tips.entrySet)
        incrementalMatrix.set(entry.key, locus, if (entry.value) 1.0 else 0.0)
      val greedy = new Greedy => [
        tipInclusionProbabilities = ReadOnlyCLMatrix::readOnly(incrementalMatrix)
      ]
      val incrementalTree = greedy.infer
      
      if (iteration % thinningPeriod === 0)
        treesWriter.write(
          "iteration" -> iteration,
          "locus" -> locus,
          "fpr" -> stat(locus).fpRate, 
          "fnr" -> stat(locus).fnRate, 
          "value" -> incrementalTree)
      iteration++
    }
  }
  static def void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}