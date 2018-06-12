package corrupt

import blang.inits.experiments.Experiment
import blang.inits.Arg
import java.io.File
import corrupt.post.CLMatrixUtils

class Distances extends Experiment {
  
  @Arg File phylo1
  @Arg File phylo2
  
  override run() {
    val mtx1 = CLMatrixUtils::fromPhylo(PerfectPhylo::parseNewick(phylo1))
    val mtx2 = CLMatrixUtils::fromPhylo(PerfectPhylo::parseNewick(phylo2))
    println("distance = " + CLMatrixUtils::distance(mtx1.matrix, mtx2.matrix)) 
  }
  
  static def main(String [] args) {
    Experiment.start(args)
  }
}