package corrupt.pre

import briefj.BriefIO
import au.com.bytecode.opencsv.CSVParser
import blang.inits.experiments.Experiment
import blang.inits.Arg
import java.io.File

class Tidify extends Experiment {
  @Arg File input
  
  override run() {
    val out = results.getAutoClosedBufferedWriter("tidy.csv")
    out.append("cells,loci,tipInclusionProbabilities\n")
    for (indexed : BriefIO::readLines(input).indexCSV(new CSVParser(" "))) {
      val locus = indexed.get("locus")
      for (entry : indexed.entrySet)
        if (entry.key != "locus")
          out.append('''«entry.key»,«locus»,«fix(entry.value)»''' + "\n")
    }
    out.close
  }
  
  public static def void main(String [] args) {
    Experiment::startAutoExit(args)
  }
  
  def static fix(String string) {
    return if (string == "0.95") 1 else 0
  }
}