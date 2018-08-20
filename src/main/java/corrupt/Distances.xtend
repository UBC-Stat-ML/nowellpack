package corrupt

import blang.inits.experiments.Experiment
import blang.inits.Arg
import java.io.File
import static corrupt.post.CLMatrixUtils.*

class Distances extends Experiment {
  
  @Arg File reference
  @Arg File guess
  
  override run() {
    val refPhylo   = PerfectPhylo::parseNewick(reference)
    val guessPhylo = PerfectPhylo::parseNewick(guess)
    results.getTabularWriter("distances") => [
      write("metric" -> "distance",          "value" -> distance(refPhylo, guessPhylo))
      write("metric" -> "locusEdgeDistance", "value" -> locusEdgeDistance(refPhylo, guessPhylo))
    ]
  }
  
  public static def void main(String [] args) {
    Experiment.start(args)
  }
}