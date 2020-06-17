package corrupt.pre

import blang.inits.experiments.Experiment
import blang.inits.Arg
import java.io.File
import corrupt.post.CLMatrixUtils
import blang.inits.DefaultValue

import static extension xlinear.MatrixExtensions.*
import corrupt.post.SimpleCLMatrix
import java.util.List
import corrupt.Locus
import corrupt.GenomeMap
import java.util.LinkedHashMap
import corrupt.Cell
import java.util.Map
import briefj.BriefMaps
import java.util.Collections
import corrupt.GenomeMap.ParsedLocus
import java.util.LinkedHashSet
import java.util.ArrayList
import blang.inits.experiments.tabwriters.factories.CSV
import briefj.BriefIO

class Filter extends Experiment {  // -> will become combine
  @Arg List<File> inputs // typically, one exec dir for negative jumps, one for positive jumps
  
  @Arg           @DefaultValue("0.05") 
  public double lowerFraction = 0.05
  
  @Arg      @DefaultValue("0.1") 
  public double maxError = 0.1
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
  
  override run() {
    var int nLowerEntriesIgnored = 0
    var int nPositiveWithinLowerEntriesIgnored = 0
    var int nCollisionsDetected = 0
    var SimpleCLMatrix data = null
    var int nLociRemaining = 0
    
    val result = new LinkedHashMap<Locus, List<Map<Cell,Double>>>
    
    for (file : inputs) {
      val matrix = CSV::csvFile(file, "binarized")
      data = CLMatrixUtils::fromCSV(matrix)
      
      for (locusLine : BriefIO::readLines(CSV::csvFile(file, "loci-statistics")).indexCSV) {
        val estimatedError = Double::parseDouble(locusLine.get("estimatedError"))
        val locus = new Locus(locusLine.get("locus"))
        val nPos = data.slice(locus).sum
        val fraction = nPos / data.cells.size
        
        if (fraction >= lowerFraction && estimatedError <= maxError) {
          nLociRemaining++
          val list = BriefMaps.getOrPutList(result, locus)
          if (!list.empty) nCollisionsDetected++
          val column = new LinkedHashMap
          for (cell : data.cells) {
            column.put(cell, data.get(cell, locus))
          }
          list.add(column)
        } 
        // those with small number of event might help estimate the FP rate
        if (fraction < lowerFraction) {
          nLowerEntriesIgnored += data.cells.size
          nPositiveWithinLowerEntriesIgnored += nPos as int
        }
      }
    }
    println("nLociRemaining = " + nLociRemaining)
    println("nCollisionsDetected = " + nCollisionsDetected)
    println("nLowerEntriesIgnored = " + nLowerEntriesIgnored)
    println("nPositiveWithinLowerEntriesIgnored = " + nPositiveWithinLowerEntriesIgnored)
    
    // build shrunk loci first
    val List<Locus> shrunkLoci = new ArrayList<Locus>
    for (entry : result.entrySet) {
      val baseLocus = entry.key
      val list = entry.value
      val parsed = new ParsedLocus(baseLocus)
      for (var int i = 0; i < list.size; i++) {
        val subLocus = GenomeMap::locus(parsed.chrString, parsed.leftOneIndexedIncl, parsed.rightOneIndexedIncl, i)
        shrunkLoci.add(subLocus)
      }
    }
    
    val allExtendedLoci = new LinkedHashSet<Locus>(shrunkLoci)
    allExtendedLoci.addAll(data.loci)
    
    for (shrinkLoci : #[true, false]) {
      val loci = if (shrinkLoci) shrunkLoci else allExtendedLoci
      val output = new SimpleCLMatrix(data.cells, loci)
      for (locus : loci) { 
        val parsed = new ParsedLocus(locus)
        val list = result.get(locus) ?: Collections::emptyList
        for (var int i = 0; i < list.size; i++) {
          val item = list.get(i)
          val subLocus = GenomeMap::locus(parsed.chrString, parsed.leftOneIndexedIncl, parsed.rightOneIndexedIncl, i)
          for (cell : data.cells) 
            output.set(cell, subLocus, item.get(cell))
        }
      }
      val name = if (shrinkLoci) "shrunk" else "full"
      CLMatrixUtils::toCSV(output, results.getFileInResultFolder("filtered-" + name + ".csv.gz")) 
    }
  }
  
}