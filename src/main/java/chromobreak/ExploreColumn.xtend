package chromobreak

import blang.inits.Arg
import blang.inits.experiments.Experiment
import briefj.BriefIO
import java.util.ArrayList
import java.util.TreeMap
import xlinear.DenseMatrix

import static xlinear.MatrixOperations.*

class ExploreColumn extends Experiment {
  
  @Arg String prefix
  @Arg String chr
  
  def static void main(String [] args){
    Experiment::startAutoExit(args)
  }
  
  override run() {
    val data = new TreeMap<Integer, DenseMatrix>
    for (line : BriefIO.readLines(prefix + ".csv").splitCSV.skip(1).filter[get(0).replaceAll("\"","") == chr]) {
      val pos = Integer.parseInt(line.get(1))
      val correctedCounts = line.subList(3, line.size).map[Double.parseDouble(it)]
      val asVector = denseCopy(correctedCounts)
      data.put(pos, asVector)
    }
    
    val sorted = new ArrayList(data.values)
    for (var int p = 2; p < sorted.size - 2; p++) {
      val writer = results.getTabularWriter(prefix + "_chr" + chr + "_deltas").child("position", p)
      val m2 = sorted.get(p - 2)
      val m1 = sorted.get(p - 1)
      val cur = sorted.get(p)
      val p1 = sorted.get(p + 1)
      val p2 = sorted.get(p + 2)
      for (c : 0 ..< cur.nEntries)
        writer.child("cell", c) => [
          write("delta" -> 1, "difference" -> (cur.get(c) - m1.get(c)))
          write("delta" -> 2, "difference" -> (p1.get(c) - m1.get(c)))
          write("delta" -> 4, "difference" -> (p2.get(c) - m2.get(c)))
          write("delta" -> -1, "difference" -> (cur.get(c) / m1.get(c)))
          write("delta" -> -2, "difference" -> ((cur.get(c) + p1.get(c)) / (m1.get(c) + m2.get(c))))
        ]

    }
  }

}