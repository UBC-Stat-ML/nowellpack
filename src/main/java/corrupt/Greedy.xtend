package corrupt

import blang.inits.experiments.Experiment
import blang.inits.Arg
import blang.inits.DefaultValue
import java.util.List
import java.util.ArrayList
import java.util.Collections
import java.util.Comparator
import corrupt.Greedy.QueuedLocus
import bayonet.distributions.Random
import java.util.Optional
import briefj.BriefIO
import corrupt.post.CellLocusMatrix
import blang.inits.experiments.tabwriters.TabularWriter

class Greedy extends Experiment {
  
  @Arg 
  public CellLocusMatrix tipInclusionProbabilities
  
  @Arg          @DefaultValue("20")
  public int reshufflePeriod = 20
  
  @Arg         @DefaultValue("20")
  public int thinningPeriod = 20
  
  @Arg 
  public Optional<Random> randomizedOrder = Optional.empty
  
  @Arg     @DefaultValue("true") // Different defaults when ran from cmd line 
  public boolean output = false  // vs called programmatically
  
  def CorruptPhylo infer() {
    val CorruptPhylo phylo = new CorruptPhylo(tipInclusionProbabilities)
    for (locus : tipInclusionProbabilities.loci)
      phylo.reconstruction.tree.collapseEdge(locus)
    val List<QueuedLocus> queue = sortQueue(null, phylo)
    var iteration = 0
    while (!queue.empty) {
      val popped = queue.pop
      val logPr = SplitSampler::maximizeInPlace(phylo.reconstruction.tree, popped,  phylo.inclusionLogProbabilities(1.0, popped))
      
      if (output)
        println("Processing locus " + (iteration+1) + "/" + tipInclusionProbabilities.loci.size + " logPr=" + logPr)
      if (iteration % thinningPeriod === 0 || queue.empty)
        writePhylo(phylo, queue, iteration, logPr, popped)      
      
      iteration++ // do before line below, already sorted at beginning
      if (iteration % reshufflePeriod === 0)
        sortQueue(queue, phylo)
    }
     
    if (output)
      BriefIO::write(results.getFileInResultFolder("tree.newick"), phylo.toString)
    return phylo
  }
  
  override run() {
    infer 
  }
  
  TabularWriter treesWriter = null
  private def void writePhylo(CorruptPhylo phylo, List<QueuedLocus> loci, int iteration, double logPr, Locus poppedLocus) {
    if (!output) 
      return;
    if (treesWriter === null) 
      treesWriter = results.getTabularWriter("trees")
    for (locus : loci)
      phylo.reconstruction.tree.addEdge(CorruptStaticUtils::root, locus.locus) 
    
    treesWriter.write(
      "iteration" -> iteration,
      "locus" -> poppedLocus,
      "logPr" -> logPr, 
      "value" -> phylo)
    
    for (locus : loci)
      phylo.reconstruction.tree.collapseEdge(locus.locus)
  }
  
  private def sortQueue(List<QueuedLocus> _old, CorruptPhylo phylo) {
    val useRandomOrder = randomizedOrder.present
    if (useRandomOrder && _old !== null)
      return _old
    val List<QueuedLocus> result = if (_old === null) {
      new ArrayList(phylo.loci.map[new QueuedLocus(it)].toList)
    } else _old
    
    if (useRandomOrder) {
      if (output)
        println("Shuffling order")
      Collections::shuffle(result, randomizedOrder.get)
    } else {
      if (output)
        println("Sorting queue")
      for (locus : result)
        locus.recomputePriority(phylo)
      Collections::sort(result, Comparator::comparing[QueuedLocus ql | ql.priority])
    } 
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
      priority = SplitSampler::maxLogConditional(phylo.reconstruction.tree, phylo.inclusionLogProbabilities(1.0, locus)) 
    }
  }
  
  static def void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}