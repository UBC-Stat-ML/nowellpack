package corrupt.pre

import blang.inits.Arg
import blang.inits.experiments.Experiment
import corrupt.post.ReadOnlyCLMatrix

/**
 * Converts tidy matrix into the format used by SCITE, 
 * i.e. in which rows are mutations and columns are loci 
 * + separate file for list of loci name.
 * Check entries are binary at same time.
 */
class Tidy2MutationCell extends Experiment {
  @Arg
  public ReadOnlyCLMatrix input
  
  override run() {
    val matrix = results.getAutoClosedBufferedWriter("matrix.csv")
    val lociNames = results.getAutoClosedBufferedWriter("matrix.geneNames")
    for (locus : input.loci) {
      lociNames.append("" + locus + "\n")
      for (cell : input.cells) {
        val value = input.get(cell, locus)
        matrix.append(
          switch (value) {
            case 0.0 : "0"
            case 1.0 : "1" 
            default : throw new RuntimeException
          }
        )
        matrix.append(" ")
      }
      matrix.append("\n")
    }
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}