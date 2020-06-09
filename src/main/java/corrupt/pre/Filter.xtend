package corrupt.pre

import blang.inits.experiments.Experiment
import blang.inits.Arg
import java.io.File
import corrupt.post.CLMatrixUtils
import java.util.HashSet
import blang.inits.DefaultValue

import static extension xlinear.MatrixExtensions.*
import corrupt.post.SimpleCLMatrix
import java.util.List
import corrupt.Locus
import corrupt.GenomeMap

class Filter extends Experiment {
  @Arg List<File> inputs // typically, one for negative jumps, one for positive jumps
  
  @Arg           @DefaultValue("0.05") 
  public double lowerFraction = 0.05

  @Arg   @DefaultValue("1.0")
  public double upperFraction = 1.0
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
  
  override run() {
    val goodLoci = new HashSet
    var int nLowerEntriesIgnored = 0
    var int nPositiveWithinLowerEntriesIgnored = 0
    var int nCollisionsDetected = 0
    var SimpleCLMatrix data = null
    for (file : inputs) {
      data = CLMatrixUtils::fromCSV(file)
      val map = new GenomeMap(data.loci)
      for (locus : data.loci) {
        val nPos = data.slice(locus).sum
        val fraction = nPos / data.cells.size
        if (fraction >= lowerFraction && fraction <= upperFraction) {
          val Locus representative = 
            if (goodLoci.contains(locus)) {
              nCollisionsDetected++
              // sometimes, we may find a positive jump 
              // at the same place as a negative jump
              // recode one as the next event (assuming that jitter will 
              // prevent this other one from being a locus)
              map.neighbors(locus, 1).get(0)
            } else {
              locus
            }
          if (goodLoci.contains(representative))
            throw new RuntimeException
          goodLoci.add(representative)
        }
        // those with small number of event might help estimate the FP rate
        if (fraction < lowerFraction) {
          nLowerEntriesIgnored += data.cells.size
          nPositiveWithinLowerEntriesIgnored += nPos as int
        }
      }
    }
    println("nCollisionsDetected = " + nCollisionsDetected)
    println("nLowerEntriesIgnored = " + nLowerEntriesIgnored)
    println("nPositiveWithinLowerEntriesIgnored = " + nPositiveWithinLowerEntriesIgnored)
    // filter
    for (shrinkLoci : #[true, false]) {
      val loci = if (shrinkLoci) goodLoci else data.loci
      val result = new SimpleCLMatrix(data.cells, loci)
      for (locus : loci) 
        for (cell : data.cells) 
          result.set(cell, locus, data.get(cell, locus))
      val name = if (shrinkLoci) "shrunk" else "full"
      CLMatrixUtils::toCSV(result, results.getFileInResultFolder("filtered-" + name + ".csv.gz")) 
    }
  }
  
}