package corrupt.pre

import blang.inits.Arg
import blang.inits.experiments.Experiment
import corrupt.post.ReadOnlyCLMatrix
import briefj.BriefIO

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
    for (locus : input.loci) {
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
    
    {
      val file = results.getFileInResultFolder("matrix.geneNames")
      BriefIO::write(file, input.loci.join("\n"))
    }
    
    {
      val file = results.getFileInResultFolder("matrix.cellNames")
      BriefIO::write(file, input.cells.join("\n"))
    }
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}