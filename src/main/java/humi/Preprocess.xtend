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
import java.util.List

class Preprocess extends Experiment {
  @Arg List<File> inputs
  
  @Arg  @DefaultValue("M") 
  String countField = "M"
  @Arg  @DefaultValue("barcode") 
  String barcodeField = "barcode"
  
  public static String HISTOGRAM_COLUMN = "histogram"
  public static String DATASET_COLUMN = "dataset"
  
  override run() {
    for (input : inputs) {
      val dataset = input.name.replaceAll("[.].*", "")
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
        lineOut.write(
          DATASET_COLUMN -> dataset,
          HISTOGRAM_COLUMN -> toString(frequencies.get(key))
        )
      }
    }
  }
  
  def static String toString(Counter<Integer> counts) {
    return counts.entries.entrySet.map["" + key + "x" + (value as double as int)].join(";")
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}