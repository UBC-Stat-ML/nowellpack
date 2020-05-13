package chromobreak

import blang.inits.experiments.Experiment
import briefj.BriefIO
import chromobreak.CombineCellPosteriors.Archipelago
import com.google.common.collect.Range
import java.io.File
import java.util.List
import java.util.LinkedHashMap
import java.util.LinkedHashSet
import java.util.Map
import binc.Command

class CorrectedDeltaAnalysis extends Experiment {
  
  def static void main(String [] args) { 
    Experiment::startAutoExit(args)
  }
  
  override run() {
  
    // create set of positions of interest for each chromo
    val pairsOfInterests = new LinkedHashSet<Pair<String,Pair<Integer,Integer>>>
  
    // load islands to know where to focus
    for (line : BriefIO::readLines("/Users/bouchard/w/nowellpack/results/all/2020-03-25-11-10-10-RVQOO6jc.exec/archipelagoes.csv").indexCSV) {
      val keep = Boolean.parseBoolean(line.get("keep"))
      if (keep) {
        val left = Integer.parseInt(line.get("leftBound"))
        val right = Integer.parseInt(line.get("rightBound"))
        val chr = line.get("chr")
        val range = Range.closedOpen(left, right)
        for (pos : CombineCellPosteriors::points(range)) {
          pairsOfInterests.add(chr -> (pos -> (pos+2))) 
//          if (pos - 1 >= 0)
//            pairsOfInterests.add(chr -> (pos -> (pos-1))) 
        }
      }
    }
    
    var Map<Pair<String,Integer>,Double> logGCs = null
    
    
    for (exec : CombineCellPosteriors::execDirs(new File("/Users/bouchard/experiments/chromobreak/links/2020-03-19-22-00-27-PmbTcHYy.exec/run"), "results/latest")) {
      val execArgs = exec.arguments
      val cellId = CombineCellPosteriors::cellId(exec)
      
      // load dataset, gc
      val dataFile = new File(execArgs.get("model.data.readCounts.dataSource"))
      if (logGCs === null) {
        val gcFile = new File(execArgs.get("model.data.source"))
        logGCs = new LinkedHashMap
        for (line : BriefIO::readLines(gcFile).indexCSV) {
          val chr = line.get("chromosomes")
          val pos = Integer.parseInt(line.get("positions"))
          val value = Double.parseDouble(line.get("value"))
          logGCs.put(chr -> pos, Math.log(value))
        }
      }
    
      // extract f1, f2 (assume it's unimodal for now)
      val f1 = univariateMean(new File(exec.execDir, "summaries/f1-summary.csv"))
      val f2 = univariateMean(new File(exec.execDir, "summaries/f2-summary.csv"))
      
      val logReads = new LinkedHashMap<Pair<String,Integer>,Double>
      for (line : BriefIO::readLines(dataFile).indexCSV) {
        val chr = line.get("chromosomes")
        val pos = Integer.parseInt(line.get("positions"))
        val value = Double.parseDouble(line.get("value"))
        logReads.put(chr -> pos, Math.log(value))        
      }
      
      for (pair : pairsOfInterests) {
        
        val chr = pair.key
        val pos1 = pair.value.key
        val pos2 = pair.value.value
        
        val key1 = chr -> pos1
        val key2 = chr -> pos2
        
        if (logReads.containsKey(key1) && logReads.containsKey(key2)) {
          val corrected1 = correct(logGCs, logReads, key1, f1, f2)
          val corrected2 = correct(logGCs, logReads, key2, f1, f2)
          val delta = corrected2 - corrected1
          val writer = results.getTabularWriter("deltas").child("cell", cellId).child("pos", pos1).child("chr", chr)
          writer.write(
            "corrected" -> true,
            "delta" -> delta
          )
          writer.write(
            "corrected" -> false,
            "delta" -> (logReads.get(key2) - logReads.get(key1))
          )
        }
        
      }
    }
    results.flushAll    
    // create plot (s?) as in EDA
    val script = '''
    require("ggplot2")
    
    data <- read.csv("«results.getFileInResultFolder("deltas.csv").absolutePath»")
    p <- ggplot(data, aes(x = delta, colour = corrected )) + 
      geom_density() + 
      facet_grid(chr + pos ~ ., scales = "free_y") + 
      theme_bw()
    ggsave("«results.getFileInResultFolder("plot.pdf").absolutePath»", p, width = 8, height = 2000, limitsize = FALSE)
    '''
    val scriptFile = results.getFileInResultFolder(".script.r")
    BriefIO.write(scriptFile, script)
    val r = Command.byName("Rscript")
    Command.call(r.appendArg(scriptFile.getAbsolutePath()))
  }
  
  def correct(Map<Pair<String, Integer>, Double> logGCs, LinkedHashMap<Pair<String, Integer>, Double> logReads, Pair<String, Integer> key, double f1, double f2) {
    val a = f2 / 2.0
    val b = f1 - 2 * a * ReadCountModel::x0
    val logRead = logReads.get(key)
    val logGC = logGCs.get(key)
    return logRead - a*logGC*logGC - b*logGC // -c ignored as we are looking at log ratios
  }
  
  def static double univariateMean(File f) {
    val list = BriefIO::readLines(f).indexCSV.map[Double.parseDouble(get("mean"))].toList
    if (list.size !== 1) throw new RuntimeException
    list.get(0)
  }
  
}