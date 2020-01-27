package chromobreak

import blang.inits.experiments.Experiment
import java.io.File
import blang.inits.Arg
import java.util.LinkedHashMap
import briefj.BriefIO

class FixGC extends Experiment {
  @Arg File csv
  override run() {
    val map = new LinkedHashMap<Pair<String,String>,String>
    for (it : BriefIO::readLines(csv).indexCSV) {
      val key = get("chromosomes") -> get("positions")
      val value = get("value")
      if (map.containsKey(key) && map.get(key) != value)
        println(value + " vs " +  map.get(key))
      map.put(key, value)
    }
    val writer = results.getTabularWriter("simplified")
    for (entry : map.entrySet)
      writer.write("chromosomes" -> entry.key.key, "positions" -> entry.key.value, "value" -> entry.value)
  }
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}