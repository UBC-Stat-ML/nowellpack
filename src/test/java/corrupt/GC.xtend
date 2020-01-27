package corrupt

import briefj.BriefIO
import briefj.BriefStrings
import blang.inits.experiments.Experiment
import blang.inits.Arg

class GC extends Experiment {
  var gcCounter = new GCStat
  var String curChr = null
  var int curBin = 0
  
  static class GCStat {
    var long total = 0
    var long gc = 0
  }
  
  @Arg(description = "Human genome in FASTA from https://www.ncbi.nlm.nih.gov/genome/guide/human/ (e.g. GRCh38)") 
  String path 
  @Arg int binSize 
  
  override run() {
    
    var long genomeLength = 0
    
    for (line : BriefIO::readLines(path)) {
      if (line.startsWith(">")) {
        dumpGCBinStats
        curChr = BriefStrings::firstGroupFromFirstMatch(".*Homo sapiens chromosome ([0-9XY]+), GRCh38.p13 Primary Assembly", line)
        if (curChr !== null) {
          println("Starting to process chromosome " + curChr)
          curBin = 0
        }
      } else if (curChr !== null) {
        val removed = line.toUpperCase.replaceAll("[ACTGN]", "")
        if (removed.length > 0) println(removed)
        val processed = line.toUpperCase.replaceAll("[^ACGTN]", "")
        genomeLength += processed.length
        gcCounter.total += processed.length
        gcCounter.gc += nGC(processed)
        if (gcCounter.total > binSize)
          dumpGCBinStats
      } else {} // skip
    }
    dumpGCBinStats
    println("Genome len = " + genomeLength)
  }
  
  def dumpGCBinStats() {
    if (gcCounter.total > 0) {
      val ratio = (gcCounter.gc as double / gcCounter.total as double)
      results.getTabularWriter("gc").write(
        "chr" -> curChr,
        "bin" -> curBin,
        "gcMean" -> ratio
      )
      curBin++
    }
    gcCounter = new GCStat
  }
  
  def static int nGC(String s) {
    return s.replaceAll("[ATN]", "").length 

  }
  
  def static void main(String [] args) {
    Experiment.startAutoExit(args) 
  }
}