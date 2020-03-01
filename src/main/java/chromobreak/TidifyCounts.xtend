package chromobreak

import briefj.BriefIO
import java.io.File
import blang.inits.experiments.Experiment
import java.util.LinkedHashMap
import java.util.TreeMap
import blang.inits.Arg
import blang.inits.DefaultValue

class TidifyCounts extends Experiment {
  
  @Arg
  public File countFile
  
  @Arg @DefaultValue("INF")
  public int nCells = Integer.MAX_VALUE
  
  override run() {
    // make sure positions are sorted
    // cell -> chr -> pos -> norm_counts
    val data = new LinkedHashMap<Pair<Integer, String>, TreeMap<Integer, Double>>()
    for (line : BriefIO.readLines(countFile).splitCSV.skip(1)) {
      val chr = line.get(0).replaceAll("\"", "")
      val pos = Integer.parseInt(line.get(1))
      for (c : 3 ..< line.size) {
        val key = Pair.of(c-3, chr)
        if (!data.containsKey(key))
          data.put(key, new TreeMap)
        val treeMap = data.get(key)
        treeMap.put(pos, Double.parseDouble(line.get(c)))
      }
    }
    
    for (key : data.keySet) {
      var int posIndex = 0
      for (pos : data.get(key).keySet) {
        val cell = key.key
        val chr = key.value
        val normCounts = data.get(key).get(pos)
        if (cell < nCells)
          results.getTabularWriter("tidy").write(
            "cell" -> cell,
            "chromosomes" -> chr,
            "positions" -> posIndex,
            "value" -> normCounts
          )
        posIndex++
      }
    }
    
  }
  
  def static void main(String [] args){
    Experiment.startAutoExit(args)
  }

}