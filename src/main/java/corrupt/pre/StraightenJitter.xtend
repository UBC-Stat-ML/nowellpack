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
import static xlinear.MatrixOperations.*
import corrupt.GenomeMap
import bayonet.math.NumericalUtils
import binc.Command
import blang.inits.experiments.tabwriters.factories.CSV
import corrupt.viz.PerfectPhyloViz
import corrupt.PerfectPhylo
import corrupt.post.ReadOnlyCLMatrix
import viz.core.Viz

class StraightenJitter extends Experiment {
  
  @Arg @DefaultValue("2") Integer neighborhoodSize = 2
  
  @Arg File input
  
  @Arg           @DefaultValue("0.05") 
  public double lowerFraction = 0.05

  @Arg           @DefaultValue("1.0")
  public double upperFraction = 1.0
  
  @Arg                @DefaultValue("0.5")
  public double posteriorThreshold = 0.5

  @Arg                @DefaultValue("Rscript")
  public Command r = Command.byName("Rscript")
  
  override run() {
    // initialize with input, will move mass around
    val result = CLMatrixUtils::fromCSV(input)
    val check = result.matrix.sum
    
    val GenomeMap lociMap = new GenomeMap(result.loci) 
    
    // order loci by prevalence
    val List<Locus> orderedLoci = orderLoci(result)
    val Set<Locus> consumed = new LinkedHashSet
    
    val goodLoci = new LinkedHashSet
    
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
        println(fraction)
        if (fraction >= lowerFraction) {
          // report merged locus summary for each above the threshold
          report(locus, result.slice(locus))
          if (fraction <= upperFraction)
            goodLoci.add(locus)
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
    
    //CLMatrixUtils::toCSV(result, results.getFileInResultFolder("binarized.csv.gz")) 
    val phylo = PerfectPhylo::parseNewick(new File("/Users/bouchard/experiments/corrupt-nextflow/results/all/2020-06-04-08-23-49-srdge3qP.exec/consensus.newick"))
    PerfectPhyloViz::visualizePerChromosome(results.getFileInResultFolder("output"), phylo, #[ReadOnlyCLMatrix::readOnly(result)], Viz::fixHeight(300))
    
    
    // heuristic: use the posterior mean of the tail of the last column (least dense ones) 
    // to set a binarization threshold
    // TODO
      
    // TODO: print the full excess probab matrix
    
    
    // binarise and print full binarised
    
    // restrict and print final result as output.csv.gz
          
    //if (viz) viz(data.matrix.copy, "final.pdf") 
    //CLMatrixUtils::toCSV(data, results.getFileInResultFolder("output.csv")) 
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