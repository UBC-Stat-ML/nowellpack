package corrupt

import org.junit.Test
import bayonet.distributions.Random
import blang.inits.experiments.Experiment
import corrupt.post.CLMatrixUtils
import corrupt.post.SimpleCLMatrix
import corrupt.post.ReadOnlyCLMatrix
import corrupt.post.CellLocusMatrix

class TestGreedy {
  
  @Test
  def void test() {
    val nCells = 1000
    val nLoci = 50
    val rand = new Random(1)
    val phylo = PerfectPhylo::generateUniform(nCells, nLoci, rand)
    
//    val matrix = new ReadOnlyCLMatrix(CLMatrixUtils::fromPhylo(phylo))  
    var _matrix = CLMatrixUtils::syntheticInclusionPrs(rand, phylo, 0.4)
    val matrix = ReadOnlyCLMatrix.readOnly(entrywiseMAP(_matrix))
    val greedy = new Greedy => [
      tipInclusionProbabilities = matrix
    ]
    val inferred = greedy.infer
    println(CLMatrixUtils::distance(phylo, inferred.reconstruction))
  }
  
  def SimpleCLMatrix entrywiseMAP(ReadOnlyCLMatrix matrix) {
    val result = new SimpleCLMatrix(matrix.cells, matrix.loci)
    for (cell : matrix.cells)
      for (locus : matrix.loci)
        result.setTip(cell, locus, round(matrix.getTipAsDouble(cell, locus)))
    return result
  }
  
  def double round(double d) {
    if (d < 0.5) return 0
    else return 1
  }
 
//  @Test
//  def void test() {
//    val nCells = 100
//    val nLoci = 50
//    val rand = new Random(1)
//    val phylo = PerfectPhylo::generateUniform(nCells, nLoci, rand)
//    
//    val matrix = new ReadOnlyCLMatrix(CLMatrixUtils::fromPhylo(phylo))  
//    val greedy = new Greedy => [
//      tipInclusionProbabilities = matrix
//    ]
//    val inferred = greedy.infer
//    println(CLMatrixUtils::distance(phylo, inferred.reconstruction))
//  }
}