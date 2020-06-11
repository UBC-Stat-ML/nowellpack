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
  
  @Arg @DefaultValue("4") Integer neighborhoodSize = 4
  
  @Arg File input
  
  @Arg                @DefaultValue("0.5")
  public double posteriorThreshold = 0.5

  @Arg                @DefaultValue("Rscript")
  public Command r = Command.byName("Rscript")
    
  @Arg                      @DefaultValue("0.05") 
  public double printDiagnosticThreshold = 0.05
  
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
    
    for (locus : orderedLoci)
      if (!consumed.contains(locus)) {
        // get left and right adjacents, (1, 2, or more of them)
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
        if (fraction >= printDiagnosticThreshold) {
          // report merged locus summary for each above the threshold
          report(locus, result.slice(locus))
        }
      }
      
    // verify we just moved mass around
    NumericalUtils::checkIsClose(result.matrix.sum, check)
    
    for (locus : result.loci) {
      for (cell : result.cells) {
        val binarized = if (result.get(cell, locus) >= posteriorThreshold) 1.0 else 0.0
        result.set(cell, locus, binarized)
      }
    }
    CLMatrixUtils::toCSV(result, results.getFileInResultFolder("binarized.csv.gz")) 
  }
  
  def report(Locus locus, Matrix vector) {
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
      
      p <- ggplot(data, aes(x = expectedCount)) + geom_density() + xlim(0, 3) + geom_rug() + theme_bw()
      ggsave("«locusSpecificResults.getFileInResultFolder("" + locus + ".pdf")»", p)
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