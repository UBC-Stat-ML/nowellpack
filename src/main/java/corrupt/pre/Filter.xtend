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
    var int nLowerEntriesIgnored = 0
    var int nPositiveWithinLowerEntriesIgnored = 0
    var int nCollisionsDetected = 0
    var SimpleCLMatrix data = null
    
    val result = new LinkedHashMap<Locus, List<Map<Cell,Double>>>
    
    for (file : inputs) {
      data = CLMatrixUtils::fromCSV(file)
      for (locus : data.loci) {
        val nPos = data.slice(locus).sum
        val fraction = nPos / data.cells.size
        
        if (fraction >= lowerFraction && fraction <= upperFraction) {
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