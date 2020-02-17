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
    output.println("map_key_1,value,sample") // log gc , log count, then need to add dummy sample column to make FitPlot r code work
    for (i : 0 ..< logReads.size)
      if (!Double.isNaN(logGCs.get(i)) && !Double.isNaN(logReads.get(i)) )
        output.println("" + logGCs.get(i) + "," + logReads.get(i) + ",0")
    output.close
    
    createPlot(
      new FitPlot(readCountModelCsvFile, types, this, rawData),
      blangExecutionDirectory.get
    )
    
    // fit histogram
    fitHistogram(rawData)
  }
  
  def fitHistogram(File rawData) {
    
    val verticalLines = '''geom_vline(xintercept = c(«(1..SingleCellHMM::nStates).map[Math::log(it)].join(",")»))'''
    
    val f0s = getFs(0)
    val f1s = getFs(1)
    val f2s = getFs(2)
    val nSamples = f0s.size
    val thin = nSamples / 10 // create 10 plots to show spread
    for (i : 0 ..< nSamples)
      if (i % thin == 0) {
        val f0 = f0s.get(i)
        val f1 = f1s.get(i)
        val f2 = f2s.get(i)
        
        val a = f2 / 2.0
        val b = f1 - 2 * a * ReadCountModel::x0
        val c = f0 - a * ReadCountModel::x0 * ReadCountModel::x0 - b * ReadCountModel::x0
        
        val output = results.getFileInResultFolder("fit-hist-" + i + "." + imageFormat)
        
        val script = '''
        require("ggplot2")
        require("dplyr")
        
        data <- read.csv("«rawData»")
        names(data)[names(data) == 'map_key_1'] <- 'logGC'
        names(data)[names(data) == 'value'] <- 'logReads'
        
        data <- data %>% mutate(transformed = logReads - «a»*logGC^2 - «b»*logGC - «c»)
        
        p <- ggplot(data, aes(x = transformed)) + geom_histogram(bins = 100)  + «verticalLines» + theme_bw()
        ggsave(filename = "«output.absolutePath»", plot = p, width = 10, height = 4, limitsize = FALSE) 
        '''
        
        callR(results.getFileInResultFolder(".fit-hist-script-" + i + ".r"), script)
      }
    
  }
  
  def List<Double> getFs(int i) {
    val sampleDir = new File(blangExecutionDirectory.get, Runner::SAMPLES_FOLDER)
    val file = new File(sampleDir, "f" + i + ".csv")
    BriefIO::readLines(file).indexCSV.map[Double.parseDouble(get("value"))].toList
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