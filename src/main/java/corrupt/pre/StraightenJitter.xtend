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
import viz.components.MatrixViz
import viz.core.Viz
import xlinear.Matrix

import static extension xlinear.MatrixExtensions.*
import corrupt.GenomeMap

class StraightenJitter extends Experiment {
  @Arg @DefaultValue("2") Integer neighborhoodSize = 2
  @Arg File input
  @Arg  @DefaultValue("false")
  public boolean viz = false
  
  override run() {
    val data = CLMatrixUtils::fromCSV(input)
    val GenomeMap lociMap = new GenomeMap(data.loci) 
    
    // viz initial matrix
    if (viz) viz(data.matrix.copy, "original.pdf")  
    
    // order loci by prevalence
    val List<Locus> orderedLoci = orderLoci(data)
    val Set<Locus> consumed = new LinkedHashSet
    
    for (locus : orderedLoci)
      if (!consumed.contains(locus)) {
        // get left and right adjacents, (1 or 2 of them)
        val List<Locus> neighbors = lociMap.neighbors(locus, neighborhoodSize)
        for (neighbor : neighbors) 
          if (!consumed.contains(neighbor)) {
            for (cell : data.cells) {
              val neighborValue = data.get(cell, neighbor)
              data.set(cell, neighbor, 0.0)
              if (neighborValue === 1.0) {
                data.set(cell, locus, 1.0) // perform logical or
              }
            }
            consumed.add(neighbor)
          }
      }
    if (viz) viz(data.matrix.copy, "final.pdf") 
    CLMatrixUtils::toCSV(data, results.getFileInResultFolder("output.csv")) 
  }
  
  def static orderLoci(SimpleCLMatrix data) {
    val result = new ArrayList(data.loci)
    Collections::sort(result, Comparator::comparing[Locus locus | data.slice(locus).sum].reversed)
    return result
  }
  
  def viz(Matrix m, String name) {
    new MatrixViz(m, MatrixViz::greyScale, Viz::fixWidth(500)).output(results.getFileInResultFolder(name))
  }
  
  def static void main(String [] args) {
    Experiment::start(args)
  }
  
}