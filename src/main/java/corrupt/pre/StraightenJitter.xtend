package corrupt.pre

import blang.inits.Arg
import blang.inits.DefaultValue
import blang.inits.experiments.Experiment
import corrupt.Locus
import corrupt.post.CLMatrixUtils
import corrupt.post.SimpleCLMatrix
import java.io.File
import java.util.ArrayList
import java.util.Collections
import java.util.Comparator
import java.util.LinkedHashSet
import java.util.List
import java.util.Set
import xlinear.Matrix
import static extension xbinc.Extensions.*


import static extension xlinear.MatrixExtensions.*
import corrupt.GenomeMap
import bayonet.math.NumericalUtils
import binc.Command
import blang.inits.experiments.tabwriters.factories.CSV

class StraightenJitter extends Experiment {
  
  @Arg @DefaultValue("4") int maxNeighborhoodSize = 4
  
  @Arg @DefaultValue("0.01") double neighbourIncreaseThreshold = 0.01
  
  @Arg @DefaultValue("1.1") double maxPerCellPosterior = 1.1
  
  @Arg File input
  
  @Arg                      @DefaultValue("0.05") 
  public double printDiagnosticThreshold = 0.05

  @Arg                @DefaultValue("Rscript")
  public Command r = Command.byName("Rscript")
  
  override run() {
    // initialize with input, will move mass around
    val result = CLMatrixUtils::fromCSV(input)
    val check = result.matrix.sum
    
    val GenomeMap lociMap = new GenomeMap(result.loci) 
    if (!lociMap.lociAdjacent)
      throw new RuntimeException("We assume all loci in a chromosome are adjacent after proper sorting. Are you running on a subset of the loci?")
    
    // order loci by prevalence
    val List<Locus> orderedLoci = orderLoci(result)
    val Set<Locus> consumed = new LinkedHashSet
    
    
    var id = 0
    for (locus : orderedLoci)
      if (!consumed.contains(locus)) {
        val neighborhoodSize = neighbouroodSize(result, lociMap, locus, consumed)
        val List<Locus> neighbors = lociMap.neighbors(locus, neighborhoodSize)
        var sumWithNeighbors = result.slice(locus).sum
        for (neighbor : neighbors) 
          if (!consumed.contains(neighbor)) {
            sumWithNeighbors += result.slice(neighbor).sum
            for (cell : result.cells) {
              val neighborValue = result.get(cell, neighbor)
              result.set(cell, neighbor, 0.0)
              result.increment(cell, locus, neighborValue)
            }
            consumed.add(neighbor)
          }
        val fraction = sumWithNeighbors / result.cells.size
        if (fraction > 1.0) throw new RuntimeException
        if (fraction >= printDiagnosticThreshold) {
          // report merged locus summary for each above the threshold
          report(id, locus, result.slice(locus))
        }
        val estimatedError = estimatedError(result.slice(locus))
        val estimatedPrevalence = estimatedPrev(result.slice(locus))
        results.getTabularWriter("loci-statistics").write(
          "id" -> id,
          "locus" -> locus,
          "estimatedError" -> estimatedError,
          "estimatedPrevalence" -> estimatedPrevalence,
          "neighborhoodSize" -> neighborhoodSize
        )
        id++
      }
      
    NumericalUtils::checkIsClose(result.matrix.sum, check)
    
    for (locus : result.loci) {
      for (cell : result.cells) {
        val binarized = if (result.get(cell, locus) >= 0.5) 1.0 else 0.0
        result.set(cell, locus, binarized)
      }
    }
    CLMatrixUtils::toCSV(result, results.getFileInResultFolder("binarized.csv.gz")) 
  }
  
  def estimatedError(Matrix matrix) {
    var sum = 0.0
    for (i : 0 ..< matrix.nEntries) {
      val entry = matrix.get(i)
      sum += Math.min(entry, 1.0 - entry)
    }
    return sum / matrix.nEntries
  }
  
  def estimatedPrev(Matrix matrix) {
    var sum = 0.0
    for (i : 0 ..< matrix.nEntries) {
      val entry = matrix.get(i)
      sum += if (entry > 0.5) 1 else 0
    }
    return sum / matrix.nEntries
  }
  
  /**
   * Increase while sum of delta posterior increases, but stopping if either of these occur:
   * - a cell is detected in which more than one event occurs
   * - we hit a consumed locus
   */
  def int neighbouroodSize(SimpleCLMatrix matrix, GenomeMap map, Locus locus, Set<Locus> consumed) {
    var lastFraction = -1.0
    for (i : 0 .. maxNeighborhoodSize) {
      val currentSum = sumWithNeighbors(matrix, map, locus, consumed, i)
      if (currentSum === null) {
        if (i == 0) throw new RuntimeException
        return i-1
      }
      val currentFraction = currentSum as double / matrix.cells.size
      val improvement = currentFraction - lastFraction
      if (improvement < 0.0) throw new RuntimeException
      if (improvement < neighbourIncreaseThreshold) return i-1
    }
    return maxNeighborhoodSize
  }
  
  // return null if detected two events or neighbours intersects with 
  def Double sumWithNeighbors(SimpleCLMatrix matrix, GenomeMap map, Locus locus, Set<Locus> consumed, int size) {
    val neighbours = map.neighbors(locus, size)
    if (!Collections.disjoint(neighbours, consumed))
      return null
    var sum = 0.0
    for (cell : matrix.cells) {
      var current = matrix.get(cell, locus)
      for (neighbour : neighbours) 
        current += matrix.get(cell, neighbour)
      if (current > maxPerCellPosterior) return null
      sum += current
    }
    return sum
  }
  
  def report(int code, Locus locus, Matrix vector) {
    val locusSpecificResults = results.child("locus-reports")
    val writer = locusSpecificResults.getTabularWriter("." + locus.toString)
    for (i : 0 ..< vector.nEntries) {
      writer.write(
        "cellIndex" -> i,
        "expectedCount" -> vector.get(i)
      )
    }
    writer.close 
    val script = results / ("." + locus + ".r")  < '''
      require("ggplot2")
      data <- read.csv("«CSV::csvFile(locusSpecificResults.resultsFolder, "." + locus.toString)»")
      
      p <- ggplot(data, aes(x = expectedCount)) + geom_density() + xlim(0, 1) + geom_rug() + theme_bw()
      ggsave("«locusSpecificResults.getFileInResultFolder("id" + code + "_" + locus + ".pdf")»", p)
    ''' 
    !(r + script)
  }
  
  def static orderLoci(SimpleCLMatrix data) {
    val result = new ArrayList(data.loci)
    Collections::sort(result, Comparator::comparing[Locus locus | data.slice(locus).sum].reversed)
    return result
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}