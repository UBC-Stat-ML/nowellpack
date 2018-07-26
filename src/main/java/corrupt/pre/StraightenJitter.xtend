package corrupt.pre

import blang.inits.Arg
import blang.inits.DefaultValue
import blang.inits.experiments.Experiment
import briefj.BriefMaps
import corrupt.Locus
import corrupt.post.CLMatrixUtils
import corrupt.post.SimpleCLMatrix
import java.io.File
import java.util.ArrayList
import java.util.Collections
import java.util.Comparator
import java.util.LinkedHashMap
import java.util.LinkedHashSet
import java.util.List
import java.util.Map
import java.util.Set
import java.util.TreeSet
import viz.components.MatrixViz
import viz.core.Viz
import xlinear.Matrix

import static extension xlinear.MatrixExtensions.*
import org.jgraph.graph.CellMapper

class StraightenJitter extends Experiment {
  @Arg @DefaultValue("2") Integer neighborhoodSize = 2
  @Arg File input
  
  override run() {
    val data = CLMatrixUtils::fromCSV(input)
    val Map<String,TreeSet<Locus>> lociMap = lociMap(data.loci) // chr -> map
    
    // viz initial matrix
    viz(data.matrix, "original.pdf") 
    
    // order loci by prevalence
    val List<Locus> orderedLoci = orderLoci(data)
    val Set<Locus> consumed = new LinkedHashSet
    
    for (locus : orderedLoci)
      if (!consumed.contains(locus)) {
        // get left and right adjacents, (1 or 2 of them)
        val List<Locus> neighbors = neighbors(lociMap, locus)
        for (neighbor : neighbors) 
          if (!consumed.contains(neighbor)) {
            for (cell : data.cells) {
              val neighborValue = data.getTipAsDouble(cell, neighbor)
              data.setTip(cell, neighbor, 0.0)
              if (neighborValue === 1.0) {
                data.setTip(cell, locus, 1.0) // or
              }
            }
            consumed.add(neighbor)
          }
      }
    viz(data.matrix, "final.pdf") 
    CLMatrixUtils::toCSV(data, results.getFileInResultFolder("output.csv")) 
  }
  
  def lociMap(Set<Locus> loci) {
    val Map<String,TreeSet<Locus>> result = new LinkedHashMap
    for (locus : loci) {
      val parsed = locus.toString.split("_")
      val chr = parsed.get(1)
      val current = BriefMaps.getOrPut(result, chr, new TreeSet<Locus>(Comparator::comparing[Locus l | Integer.parseInt(l.toString.split("_").get(2))]))
      current.add(locus)
    }
    return result
  }
  
  def neighbors(Map<String,TreeSet<Locus>> maps, Locus locus) {
    val result = new ArrayList<Locus>
    val parsed = locus.toString.split("_")
    val chr = parsed.get(1)
    val map = maps.get(chr)
    var Locus cur = null
    cur = locus; for (i : 1..neighborhoodSize) { if (cur !== null) cur = map.higher(cur); if (cur !== null) result.add(cur) }
    cur = locus; for (i : 1..neighborhoodSize) { if (cur !== null) cur = map.lower(cur);  if (cur !== null) result.add(cur) }
    return result
  }
  
  def static orderLoci(SimpleCLMatrix data) {
    val result = new ArrayList(data.loci)
    Collections::sort(result, Comparator::comparing[Locus locus | -data.slice(locus).sum])
    return result
  }
  
  def viz(Matrix m, String name) {
    new MatrixViz(m, MatrixViz::greyScale, Viz::fixWidth(500)).output(results.getFileInResultFolder(name))
  }
  
  def static void main(String [] args) {
    Experiment::start(args)
  }
  
}