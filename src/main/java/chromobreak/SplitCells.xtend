package chromobreak

import blang.inits.experiments.Experiment
import blang.inits.Arg
import java.io.File
import briefj.BriefIO

class SplitCells extends Experiment {
  @Arg File csv
  override run() {
    for (line : BriefIO::readLines(csv).indexCSV) 
      results
        .getTabularWriter(line.get("cell"))
        .write(
          line.entrySet.filter[key != "cell"].map[Pair.of(key,value)]
        )
  }
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}