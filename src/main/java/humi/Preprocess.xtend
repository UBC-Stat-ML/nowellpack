package humi

import blang.inits.experiments.Experiment
import blang.inits.Arg
import java.io.File
import briefj.BriefIO
import java.util.Map
import briefj.collections.Counter
import static extension briefj.BriefMaps.getOrPut
import blang.inits.DefaultValue
import blang.inits.experiments.tabwriters.factories.CSV
import java.util.LinkedHashMap

class Preprocess extends Experiment {
  @Arg File input
  
  @Arg  @DefaultValue("M") 
  String countField = "M"
  @Arg  @DefaultValue("barcode") 
  String barcodeField = "barcode"
  
  override run() {
    val frequencies = new LinkedHashMap<Map<String,String>, Counter<Integer>>
    for (line : BriefIO::readLines(input).indexCSV) {
      val count = Integer.parseInt(line.remove(countField))
      if (count > 0) {
        line.remove("")
        line.remove(barcodeField)
        frequencies.getOrPut(line, new Counter).incrementCount(count, 1)
      }
    }
    
    val out = results.getTabularWriter("output", new CSV) 
    for (key : frequencies.keySet) {
      var lineOut = out
      for (entry : key.entrySet) {
        lineOut = lineOut.child(entry.key, entry.value)
      }
      lineOut.write("histogram" -> toString(frequencies.get(key)))
    }
  }
  
  def static String toString(Counter<Integer> counts) {
    return counts.entries.entrySet.map["" + key + "x" + (value as double as int)].join(";")
  }
  
  def public static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}