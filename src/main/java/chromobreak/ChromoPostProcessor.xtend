package chromobreak

import blang.runtime.internals.DefaultPostProcessor
import java.io.File
import blang.runtime.Runner
import java.util.Map
import blang.inits.experiments.tabwriters.TidySerializer
import briefj.BriefIO
import java.util.List

class ChromoPostProcessor extends DefaultPostProcessor {
  
  public static SingleCellData data = null // auto-set the above from SingleCell if ran as postprocessor
  
  override run() {    
    val outputDir = results.getFileInResultFolder("chromoplots")
    val readCountModelName = "readCountModel"
    val sampleDir = new File(blangExecutionDirectory.get, Runner::SAMPLES_FOLDER)
    val readCountModelCsvFile = new File(sampleDir, readCountModelName + ".csv") 
    val types = TidySerializer::types(readCountModelCsvFile)
    types.remove("map_key_1") // don't auto-facet this
    
    val rawData = new File(outputDir, ".raw.csv")
    val output = BriefIO::output(rawData)
    output.println("chromosomes,positions,logGC,logReads,sample") // log gc , log count, then need to add dummy sample column to make FitPlot r code work
    for (chr : data.chromosomes.indices)
      for (pos : data.positions.indices(chr)) {
        val logRead = Math::log(data.readCounts.get(chr, pos).intValue)
        val logGC = Math::log(data.gcContents.get(chr, pos).doubleValue)
        if (!Double.isNaN(logGC) && !Double.isNaN(logRead))
         output.println("" + chr.key + "," + pos.key + "," + logGC + "," + logRead + ",0")
      }
    output.close
    
    createPlot(
      new FitPlot(readCountModelCsvFile, types, this, rawData),
      outputDir
    )
    
    // state paths
    val _hmms = new File(sampleDir, "hmms.csv")
    if (_hmms.exists) {
      // workaround: if we leave this in samples, this will create huge facet plots, 
      // which are not most useful representation; we move them up to avoid this
      val hmms = new File(outputDir, _hmms.name)
      _hmms.renameTo(hmms)
      paths(hmms)
      goodnessOfFit(outputDir, rawData, hmms)
    } else 
      goodnessOfFit(outputDir, rawData, null)
    
    super.run
  }
    
  public static int nStates = 12
  
  def String means() {
    (1..nStates).map[Math::log(it)].join(",")
  }
  
  def goodnessOfFit(File outputDir, File rawData, File paths) {
    
    val verticalLines = '''geom_vline(xintercept = c(«means»))'''
    
    val f0s = getFs(0)
    val f1s = getFs(1)
    val f2s = getFs(2)
    val sds = getUnivariateSampleList("sd")
    val slopes = getUnivariateSampleList("sdSlope")
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
        
        { // histogram diagnostic
          val output = new File(outputDir, "fit-hist-" + i + "." + imageFormat)
          val script = '''
          require("ggplot2")
          require("dplyr")
          
          data <- read.csv("«rawData»")
          
          data <- data %>% mutate(transformed = logReads - «a»*logGC^2 - «b»*logGC - «c»)
          
          p <- ggplot(data, aes(x = transformed)) + geom_histogram(bins = 100)  + «verticalLines» + theme_bw()
          ggsave(filename = "«output.absolutePath»", plot = p, width = 10, height = 4, limitsize = FALSE) 
          '''
          
          callR(new File(outputDir, ".fit-hist-script-" + i + ".r"), script)
        }
        
        if (paths !== null) { // path diagnostic
          val output = new File(outputDir, "fit-path-" + i + "." + imageFormat)
          
          val script = '''
          require("ggplot2")
          require("dplyr")
          
          data <- read.csv("«rawData»")
          data <- data %>% mutate(value = exp(logReads - «a»*logGC^2 - «b»*logGC - «c»))
          
          data2 <- read.csv("«paths»") %>% filter(sample == «i»)
          names(data2)[names(data2) == 'map_key_0'] <- 'chromosomes'

          p <- ggplot(data, aes(x = positions, y = value)) +
                          geom_line(data = data2) + 
                          geom_point(data = data, alpha = 0.1) + 
                          theme_bw() +
                          facet_grid(. ~ chromosomes) + 
                          scale_y_continuous(breaks=seq(0,«nStates»,1), limits = c(0,«nStates»), minor_breaks = NULL) +
                          theme(legend.position = "none") +
                          ylab("states") +
                          xlab("positions") 
          ggsave(filename = "«output.absolutePath»", plot = p, width = 100, height = 5, limitsize = FALSE)
          '''
          callR(new File(outputDir, ".fit-gof-script-" + i + ".r"), script)
          
          val gofDir = new File(outputDir, "gof")
          for (j : 1 .. nStates) {
            val out2 = new File(gofDir, "iteration=" + i + ",state=" + j + "." + imageFormat)
            val s2 = '''
              require("ggplot2")
              require("dplyr")
              
              data <- read.csv("«rawData»")
              data <- data %>% mutate(value = exp(logReads - «a»*logGC^2 - «b»*logGC - «c»))
              
              data2 <- read.csv("«paths»") %>% filter(sample == «i», value == «j»)
              names(data2)[names(data2) == 'map_key_0'] <- 'chromosomes'              
              
              merged <- inner_join(data, data2, by = c("chromosomes", "positions")) 
              
              p <- ggplot(merged, aes(x = log(value.x))) +
                             stat_function(fun = dnorm, n = 1000, args = list(mean = log(«j»), sd = «sds.get(i) + slopes.get(i) * j»), colour = "red") + 
                             «verticalLines» + 
                             geom_density() + 
                             theme_bw() +
                             ylab("transformed read counts") +
                             xlab("density (black:data, red:model)") 
              ggsave(filename = "«out2.absolutePath»", plot = p, width = 10, height = 20, limitsize = FALSE)
            '''
            callR(new File(gofDir, ".iteration=" + i + ",state=" + j + ".r"), s2)
          }
        }
      }
  }
  
  def List<Double> getFs(int i) {
    getUnivariateSampleList("f" + i)
  }
  
  def List<Double> getUnivariateSampleList(String name) {
    val sampleDir = new File(blangExecutionDirectory.get, Runner::SAMPLES_FOLDER)
    val file = new File(sampleDir, name + ".csv")
    BriefIO::readLines(file).indexCSV.map[Double.parseDouble(get("value"))].toList
  }
  
  def paths(File posteriorSamples) {
    
    val output = results.getFileInResultFolder("paths." + imageFormat)
    val script = '''
    require("ggplot2")
    require("dplyr")
    
    data <- read.csv("«posteriorSamples»")
    n_samples <- max(data$«Runner.sampleColumn»)
    cut_off <- n_samples * «burnInFraction»
    data <- subset(data, «Runner.sampleColumn» > cut_off)
    p <- ggplot(data[!is.na(data$value),], aes(x = positions, y = value, colour = factor(«Runner::sampleColumn»))) +
                    geom_line(alpha = 0.1) + 
                    theme_bw() +
                    facet_grid(. ~ map_key_0) + 
                    scale_y_continuous(breaks=seq(0,«nStates»,1), limits = c(0,«nStates»), minor_breaks = NULL) +
                    theme(legend.position = "none") +
                    ylab("states") +
                    xlab("positions") 
    ggsave(filename = "«output.absolutePath»", plot = p, width = 100, height = 5, limitsize = FALSE) 
    '''
    
    callR(results.getFileInResultFolder(".paths.r"), script)
    
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
      names(df2)[names(df2) == 'logGC'] <- 'map_key_1'
      names(df2)[names(df2) == 'logReads'] <- 'value'
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

}