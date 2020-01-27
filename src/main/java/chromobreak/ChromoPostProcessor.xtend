package chromobreak

import blang.runtime.internals.DefaultPostProcessor
import java.io.File
import blang.runtime.Runner
import java.util.Map
import blang.inits.experiments.tabwriters.TidySerializer
import briefj.BriefIO
import java.util.List
import java.util.ArrayList

class ChromoPostProcessor extends DefaultPostProcessor {
    
  
  override run() {
    super.run
    
    val readCountModelName = "readCountModel"
    val sampleDir = new File(blangExecutionDirectory.get, Runner::SAMPLES_FOLDER)
    val readCountModelCsvFile = new File(sampleDir, readCountModelName + ".csv")
    val types = TidySerializer::types(readCountModelCsvFile)
    types.remove("map_key_1") // don't auto-facet this
    
    val rawData = new File(blangExecutionDirectory.get, ".raw.csv")
    val output = BriefIO::output(rawData)
    output.println("map_key_1,value,sample") // log gc , log count
    for (i : 0 ..< logReads.size)
      if (!Double.isNaN(logGCs.get(i)) && !Double.isNaN(logReads.get(i)) )
        output.println("" + logGCs.get(i) + "," + logReads.get(i) + ",0")
    output.close
    
    createPlot(
      new FitPlot(readCountModelCsvFile, types, this, rawData),
      blangExecutionDirectory.get
    )
  }
  
  static class FitPlot extends GgPlot {
    val File rawData
    new(File posteriorSamples, Map<String, Class<?>> types, DefaultPostProcessor processor, File rawData) {
      super(posteriorSamples, types, processor)
      this.rawData = rawData
    }
    override ggCommand() {
      return '''
      «removeBurnIn»
      df2 <- read.csv("«rawData.absolutePath»")
      p <- ggplot(data, aes(x = map_key_1, y = value, colour = factor(«Runner::sampleColumn»))) +
        geom_line(data = data, alpha = 0.1) + «facetString»
        geom_point(data = df2, alpha = 0.1) +
        theme_bw() +
        theme(legend.position = "none") +
        ylab("log read count") +
        xlab("log GC contents") 
      '''
    }
  }
  
  static List<Double> logReads = new ArrayList
  static List<Double> logGCs = new ArrayList
  def static void addToPlot(double[] _logReads, double[] _logGCs) {
    logReads.addAll(_logReads)
    logGCs.addAll(_logGCs)
  }
}