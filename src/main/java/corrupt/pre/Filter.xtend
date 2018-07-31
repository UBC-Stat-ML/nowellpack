package corrupt.pre

import blang.inits.experiments.Experiment
import blang.inits.Arg
import java.io.File
import corrupt.post.CLMatrixUtils
import java.util.HashSet
import blang.inits.DefaultValue

import static extension xlinear.MatrixExtensions.*
import corrupt.post.SimpleCLMatrix

class Filter extends Experiment {
  @Arg File input
  
  @Arg           @DefaultValue("0.0") 
  public double lowerFraction = 0.0

  @Arg   @DefaultValue("1.0")
  public double upperFraction = 1.0
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
  
  override run() {
    val data = CLMatrixUtils::fromCSV(input)
    // find set of loci
    val goodLoci = new HashSet
    for (locus : data.loci) {
      val fraction = (data.slice(locus).sum as double) / data.cells.size
      if (fraction >= lowerFraction && fraction <= upperFraction)
        goodLoci.add(locus)
    }
    // filter
    val result = new SimpleCLMatrix(data.cells, goodLoci)
    for (locus : goodLoci) 
      for (cell : data.cells) 
        result.set(cell, locus, data.get(cell, locus))
    CLMatrixUtils::toCSV(result, results.getFileInResultFolder("filtered.csv")) 
  }
  
}