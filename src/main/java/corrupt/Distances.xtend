package corrupt

import blang.inits.experiments.Experiment
import blang.inits.Arg
import java.io.File
import corrupt.post.CLMatrixUtils
import xlinear.MatrixOperations
import bayonet.distributions.Random
import java.util.stream.Collectors

class Distances extends Experiment {
  
  @Arg File reference
  @Arg File guess
  
  override run() {
    val refPhylo = PerfectPhylo::parseNewick(reference)
    val refMtx = CLMatrixUtils::fromPhylo(refPhylo)
    val guessMtx = CLMatrixUtils::fromPhylo(PerfectPhylo::parseNewick(guess))
    println("distance = " + CLMatrixUtils::distance(refMtx, guessMtx)) 
    println("allZeroBaseline = " + CLMatrixUtils::distance(refMtx.matrix, MatrixOperations::dense(refMtx.matrix.nRows, refMtx.matrix.nCols)))
    
    val summaryStats = (0..10).toList.stream.collect(Collectors.summarizingDouble[
      val randomPhylo = new PerfectPhylo(refPhylo.cells, refPhylo.loci)  
      randomPhylo.sampleUniform(new Random(1))
      CLMatrixUtils::distance(refMtx, CLMatrixUtils::fromPhylo(randomPhylo))
    ])
    println("randomBaseline = " + summaryStats.average)
  }
  
  public static def void main(String [] args) {
    Experiment.start(args)
  }
}