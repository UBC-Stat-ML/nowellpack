package humi.posets

import blang.inits.experiments.Experiment
import blang.inits.parsing.Posix
import humi.HumiStaticUtils
import java.io.File
import java.util.LinkedHashSet
import briefj.BriefIO
import humi.freq.DeltaMethod
import blang.inits.Arg
import blang.inits.DefaultValue
import java.util.LinkedHashMap

class Intervals2Poset extends Experiment {
  
  @Arg @DefaultValue("sgrna")
  public String sgRNAField = "sgrna"
  
  @Arg @DefaultValue("gene")
  public String geneAField = "gene"
  
  @Arg @DefaultValue("0.5")
  public double threshold = 0.5
  
  @Arg public File intervalsCSVFile
  
  override run() {
    val poset = loadPoset
    Posets::dotExporter(poset).export(results.getFileInResultFolder("hasse.dot"))
  }
  
  def Poset<String> loadPoset() {
    if (threshold < 0.0) throw new RuntimeException
    
    val leftIntervals = new LinkedHashMap<String,Double>
    val rightIntervals = new LinkedHashMap<String,Double>
    val names = new LinkedHashSet
    
    for (line : BriefIO::readLines(intervalsCSVFile).indexCSV) {
      val name = line.get(geneAField) + " (" + line.get(sgRNAField) + ")"
      if (names.contains(name)) throw new RuntimeException
      names.add(name)
      val left = Double.parseDouble(line.get(DeltaMethod.Columns::logRatioLeftBound.toString))
      val right = Double.parseDouble(line.get(DeltaMethod.Columns::logRatioRightBound.toString))
      
      if (right < -threshold || left > threshold) {
        leftIntervals.put(name, left)
        rightIntervals.put(name, right)
      }
    }
     
    val controls = "controls"
    leftIntervals.put(controls, 0.0)
    rightIntervals.put(controls, 0.0)
    
    val poset = new Poset<String>() {
      override compare(String first, String second) {
        if (rightIntervals.get(first) < leftIntervals.get(second)) return -1
        if (rightIntervals.get(second) < leftIntervals.get(first)) return 1
        return null
      }
      override objects() {
        return leftIntervals.keySet
      }
    }
    
    return GraphPoset.from(poset) 
  }
  
  def static void main(String [] args) {
    Experiment::start(args, Posix::parse(args), HumiStaticUtils::parsingConfigs)
  }
}