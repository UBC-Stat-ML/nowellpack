package corrupt

import blang.inits.experiments.Experiment
import blang.inits.Arg
import corrupt.post.ReadOnlyCLMatrix
import blang.inits.DefaultValue
import java.util.List
import java.util.ArrayList
import java.util.Collections
import java.util.Comparator
import briefj.BriefIO

class Greedy extends Experiment {
  
  @Arg ReadOnlyCLMatrix tipInclusionProbabilities
  
  @Arg   @DefaultValue("10")
  int reshufflePeriod = 10
  
  override run() {
    val CorruptPhylo phylo = new CorruptPhylo(tipInclusionProbabilities)
    val List<QueuedLocus> queue = sortQueue(null, phylo)
    var iteration = 0
    while (!queue.empty) {
      println("Processing locus " + (iteration+1) + "/" + tipInclusionProbabilities.loci.size)
      val popped = queue.pop
      SplitSampler::maximize(phylo.reconstruction.tree, popped,  phylo.cellInclusionLogProbabilities(1.0, popped))
      iteration++ // do before line below, already sorted at beginning
      if (iteration % reshufflePeriod == 0)
        sortQueue(queue, phylo)
    }
    val outputFile = results.getFileInResultFolder("tree.newick")
    BriefIO::write(outputFile, phylo.toString)
  }
  
  private def sortQueue(List<QueuedLocus> _old, CorruptPhylo phylo) {
    println("Sorting queue")
    val List<QueuedLocus> result = if (_old === null) {
      new ArrayList(phylo.loci.map[new QueuedLocus(it)].toList)
    } else _old
    for (locus : result)
      locus.recomputePriority(phylo)
    Collections::sort(result, Comparator::comparing[QueuedLocus ql | ql.priority]) 
    return result
  }
  
  private def pop(List<QueuedLocus> queue) {
    return queue.remove(queue.size - 1).locus
  }
  
  static class QueuedLocus {
    val Locus locus
    new (Locus l) { this.locus = l }
    var double priority = Double.NaN
    def recomputePriority(CorruptPhylo phylo) {
      priority = SplitSampler::maxLogConditional(phylo.reconstruction.tree, locus, phylo.cellInclusionLogProbabilities(1.0, locus)) 
    }
  }
  
  static def void main(String [] args) {
    Experiment::start(args)
  }
}