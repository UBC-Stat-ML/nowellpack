package chromobreak

import blang.inits.experiments.Experiment
import java.io.File
import blang.inits.Arg
import blang.xdoc.BootstrapHTMLRenderer
import blang.xdoc.components.Document
import briefj.BriefFiles
import blang.xdoc.components.Embed
import blang.inits.DefaultValue
import blang.inits.experiments.doc.ExperimentHTMLDoc.ParsedExperiment
import briefj.BriefIO
import org.eclipse.xtend.lib.annotations.Data
import java.util.LinkedHashMap
import briefj.BriefMaps
import briefj.collections.Counter
import java.util.Collections
import org.eclipse.xtext.xbase.lib.Functions.Function1
import java.util.LinkedHashSet
import com.google.common.collect.Range
import com.google.common.collect.Multimap
import com.google.common.collect.LinkedHashMultimap
import java.util.Map
import org.apache.commons.math3.stat.descriptive.SummaryStatistics
import chromobreak.datastructures.Interval
import java.util.Iterator
import chromobreak.datastructures.IntervalTree
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics
import java.util.HashMap
import java.util.List
import java.util.ArrayList
import binc.Command
import java.util.HashSet

class CombineCellPosteriors extends Experiment {
  
  @Arg File directory
  
  @Arg @DefaultValue("results/latest")
  String execInDirectory = "results/latest"
  
  @Arg(description = "Rank cells per mean island length. Remove this fraction of the long ones.")
                         @DefaultValue("0.5")
  double fractionIslandLengthCutoff = 0.5
  
  @Arg @DefaultValue("0.5")
  double burnInFraction = 0.5
  
  @Arg(description = "Mass of change points with delta >= 0 is added at each location. Only mass above that threshold will open or extend 'islands'") 
       @DefaultValue("0.02")
  double  threshold = 0.02
  
  @Arg(description = "After an island is created, sum to mass of all positions involved. Exclude if below this cut off'") 
         @DefaultValue("0.1")
  double   massCutOff = 0.1 
  
  @Arg   @DefaultValue("Rscript")
  public String rCmd = "Rscript"
  
  // ALL cell ids (some of which will be excluded; those that have intervals too wide)
  val cellIds = new LinkedHashSet<String> 
  
  // this gets populated by loadData()
  var LinkedHashMap<String,Integer> lengths = null // chr -> length
  
  val expectedNumberOfEvents = new Counter<String>
  
  @Data
  static class Key {
    val String cell
    val String chr
    val int pos
  }
  
  // create posterior summary across several cell-specific exec directories
  // returns key -> count over different deltas
  // includes ALL change points, no filtering done at this stage
  def Map<Key,Counter<Integer>> loadData() {
    // key -> counts over different deltas
    val changePoints = new LinkedHashMap<Key,Counter<Integer>>
    for (exec : execDirs) {
      val cellId = cellId(exec)
      cellIds.add(cellId)
      val paths = new File(exec.execDir, "chromoplots/hmms.csv")
      var lines = BriefIO::readLines(paths).indexCSV.toList
      if (lengths === null) {
        loadChromosomeLengths(lines)
      }
      val nIters = 1+lines.map[Integer.parseInt(get("sample"))].max
      val burnIn = burnInFraction * nIters
      val nItersUsed = nIters - burnIn
      for (i : 1 ..< lines.size - 1) {
        val state = lines.get(i).get("value")
        val pc = lines.get(i-1).get("map_key_0")
        val nc = lines.get(i+1).get("map_key_0")
        val sample = Integer.parseInt(lines.get(i).get("sample"))
        if (state == "NA" && pc == nc && sample >= burnIn) {
          val prev = Integer::parseInt(lines.get(i-1).get("value"))
          val next = Integer::parseInt(lines.get(i+1).get("value"))
          val delta = next - prev
          val key = new Key(cellId, nc, Integer.parseInt(lines.get(i).get("positions")))
          BriefMaps.getOrPut(changePoints, key, new Counter).incrementCount(delta, 1.0/nItersUsed)
          expectedNumberOfEvents.incrementCount(cellId, 1.0/nItersUsed)
        }
      }
    }
    println(expectedNumberOfEvents)
    return changePoints
  }
  
  def record(Map<Key,Counter<Integer>> changePoints) {
    for (key : changePoints.keySet) {
      val counter = changePoints.get(key)
      for (delta : counter) {
        results.getTabularWriter("changePoints").
          write(
            "cells" -> key.cell,
            "chromosomes" -> key.chr,
            "positions" -> key.pos,
            "delta" -> delta,
            "count" -> counter.getCount(delta)
          )
      }
    }
  }
  
  def loadChromosomeLengths(List<Map<String,String>> lines) {
    lengths = new LinkedHashMap<String,Integer>
    for (line : lines) {
      val chr = line.get("map_key_0")
      val len = Integer.parseInt(line.get("positions")) + 1
      val cur = lengths.get(chr) ?: 0
      lengths.put(chr, Math::max(len, cur))
    }
  }
  
  def runChecks() {
    if (fractionIslandLengthCutoff < 0.0 || fractionIslandLengthCutoff > 1.0) throw new RuntimeException
  }
  
  /* 
   * filterings done here:
   * - taking only account increases in copy numbers (or zero)
   * - remove cells where intervals are too large
   */
  def Multimap<String,Island> islands(Map<Key,Counter<Integer>> changePoints) {
    val Multimap<String,Island> islands = LinkedHashMultimap.create // cellId -> island
    val meanLenStats = new DescriptiveStatistics
    val meanLens = new HashMap<String,Double> // cellId -> mean len
    for (cellId : cellIds) {
      val islandLength = new SummaryStatistics
      for (chr : lengths.keySet) {
        var int start = -1
        var Counter<Integer> posterior = new Counter
        for (pos : 0 ..< lengths.get(chr)) {
          val key = new Key(cellId, chr, pos)
          val counter = changePoints.get(key) ?: new Counter
          val massOnCopyNumberIncrease = massOnCopyNumberIncrease(counter)
          if (massOnCopyNumberIncrease > threshold && /* This conditions ensures the last island of the chr is closed */ pos < lengths.get(chr) - 1) {
            if (start === -1) {
              // open new one
              start = pos
            } else {
              // keep extending (nothing to do)
            }
            posterior.incrementAll(counter)
          } else {
            if (start !== -1) {
              // close it
              for (delta : posterior.keySet)
                posterior.setCount(delta, posterior.getCount(delta))
              if (posterior.totalCount > massCutOff) {
                val island = new Island(cellId, chr, Range.closedOpen(start, pos), posterior)
                islandLength.addValue(island.pos.upperEndpoint - island.pos.lowerEndpoint) 
                islands.put(cellId, island)
              }
              start = -1
              posterior = new Counter
            }
          }
        }
      }
      meanLens.put(cellId, islandLength.mean)
      meanLenStats.addValue(islandLength.mean)
    }
    
    // remove cells where islands are too large
    val cutoff = meanLenStats.getPercentile(fractionIslandLengthCutoff * 100.0)
    for (cellId : cellIds) {
      val meanLen = meanLens.get(cellId)
      val keep = meanLen <= cutoff
      results.getTabularWriter("islandMeanLength").write(
        "cell" -> cellId,
        "value" -> meanLen,
        "kept" -> keep
      )
      if (!keep)
        islands.removeAll(cellId)
    }
    
    println("nCellsLeft = " + islands.nCellsLeft)
    return islands
  }
  
  def static int nCellsLeft(Multimap<String, Island> islands) {
    return islands.keySet.size
  }
  
  def Map<String,IntervalTree<Archipelago>> archipelagoes(Multimap<String,Island> islands) {
    val intervalTrees = new LinkedHashMap<String,IntervalTree<Archipelago>> // chr -> tree
    for (island : islands.values) {
      val tree = BriefMaps.getOrPut(intervalTrees, island.chr, new IntervalTree)
      val archi = new Archipelago(island.pos, island.chr, build(island))
      var newRange = island.pos 
      var newList = archi.islands     
      for (other : [tree.overlappers(archi)]) {
        tree.delete(other)
        newRange = minimumCoveringInterval(newRange, other.pos)
        newList = coalesce(newList, other.islands)
      }
      tree.insert(new Archipelago(newRange, island.chr, newList))
    } 
    
    // compute gaps between previous and next achipelogos (to estimate local intensity)
    val minGapStats = new Counter
    for (archis : intervalTrees.values) {
      val Iterable<Archipelago> iterable = [archis.iterator]
      for (archi : iterable) {
        val prev = archis.predecessor(archi)
        val next = archis.successor(archi)
        archi.gaps.add(archi.pos.lowerEndpoint - (if (prev.present) prev.get.pos.upperEndpoint else 0))
        archi.gaps.add((if (next.present) next.get.pos.lowerEndpoint else lengths.get(archi.chr)) - archi.pos.upperEndpoint)
        minGapStats.incrementCount(archi.gaps.min, 1.0)
      }
    }
    println(minGapStats) 
    return intervalTrees
  }
  
  def void record(Map<String,IntervalTree<Archipelago>> intervalTrees, Multimap<String,Island> islands, Map<Key,Counter<Integer>> changePoints) {
    val archiKeptPDF = results.child("archipelagoes")
    val archiRemPDF = results.child("archipelagoes-excluded")
    val archiCSV = results.child("archipelagoes-CSVs")
    for (archis : intervalTrees.values) {
      val Iterable<Archipelago> iterable = [archis.iterator]
      for (archi : iterable) {
        var keep = true
        var msg = new ArrayList<String>
        if (archi.islands.size === 1 || archi.islands.size === islands.nCellsLeft) {
          keep = false
          msg += "Event not informative as it covers either all cells, or only one"
        }
        results.getTabularWriter("archipelagoes").write(
          "id" -> archi.id,
          "chr" -> archi.chr,
          "leftBound" -> archi.pos.lowerEndpoint,
          "rightBound" -> archi.pos.upperEndpoint,
          "leftGrap" -> archi.gaps.get(0),
          "rightGrap" -> archi.gaps.get(1),
          "nIslands" -> archi.islands.size,
          "fractionIslands" -> archi.islands.size as double / islands.nCellsLeft,
          "minMass" -> archi.minMass,
          "maxMass" -> archi.maxMass,
          "keep" -> keep,
          "keepMessage" -> msg.join("; ")
        )
        val cells = new HashSet<String>
        
        val tabWriter = archiCSV.getTabularWriter(archi.id)
        for (cell : islands.keySet) {
          for (pos : points(archi.pos)) {
            val key = new Key(cell, archi.chr, pos) 
            val counter = changePoints.get(key)
            if (counter !== null) {
              for (delta : counter.keySet) 
                tabWriter.write(
                  "cell" -> cell,
                  "pos" -> pos,
                  "delta" -> delta,
                  "value" -> counter.getCount(delta)
                )
              cells.add(cell)  
            }
          }
        }
        archiCSV.flushAll
        val folder = (if (keep) archiKeptPDF else archiRemPDF)
        callR(archiCSV.getFileInResultFolder("." + archi.id + ".r"), '''
          require("ggplot2")
          data <- read.csv("«archiCSV.getFileInResultFolder(archi.id + ".csv").absolutePath»")
          p <- ggplot(data,aes(x=delta,y=value,fill=factor(delta))) + 
            geom_bar(stat="identity") + 
            facet_grid(cell ~ pos) + 
            theme_bw()
          ggsave("«folder.getFileInResultFolder(archi.id + ".pdf").absolutePath»", p, limitsize = FALSE, height = «3*cells.size»)
        ''')
        
        

        
      }
    }
  }
  
  override run() {
    // TODO: exclude locations with no GC contents info?
    runChecks
    val changePoints = loadData()
    record(changePoints)
    val islands = islands(changePoints)
    val archipelagoes = archipelagoes(islands)
    record(archipelagoes, islands, changePoints)
    // report
    new Website(this).renderInto(results.getFileInResultFolder("output"))
    
    /*
     * TODO:
     * - ignore events dominated by 0
     * - cf with HMM copy
     * - report total activity per cell [done]
     * - discard events involving copy number 1
     * - filtering of cell of wide intervals too aggressive??
     * 
     * - check if weird GC relations indeed got filtered out?
     * 
     * - go back to delta type approach, but this time with 
     *   starting to look like event might have "side effects"
     *   i.e. neightbors sites with higher rate. HMMs might 
     *   be problematic in this context since it might merge 
     *   groups of events
     * - hope is that removing cell-specific random effect makes peeks much more clean (right now they overlap several deltas)  
     * 
     */
  }
  
  def static points(Range<Integer> pos) {
    (pos.lowerEndpoint ..< pos.upperEndpoint)
  }
  
  def static Range<Integer> minimumCoveringInterval(Range<Integer> r1, Range<Integer> r2) {
    return Range.closedOpen(Math::min(r1.lowerEndpoint, r2.lowerEndpoint), Math::max(r1.upperEndpoint, r2.upperEndpoint))
  }
  
  def static double massOnCopyNumberIncrease(Counter<Integer> counter) {
    var sum = 0.0
    for (delta : counter.keySet.filter[it >= 0]) {
      sum += counter.getCount(delta)
    }
    return sum
  }
  
  @Data
  static class Archipelago implements Interval {
    val Range<Integer> pos
    val String chr
    val CoalescingList<Island> islands
    val List<Integer> gaps = new ArrayList
    
    def String id() { chr + "@" + pos }
    
    def double minMass() {
      islands.map[posterior.totalCount].min
    }
    
    def double maxMass() {
      islands.map[posterior.totalCount].max
    }
    
    override start() {
      pos.lowerEndpoint
    }
    
    override end() {
      pos.upperEndpoint
    }
  }
  
  def static <T> CoalescingList<T> coalesce(CoalescingList<T> l1, CoalescingList<T> l2) {
    l1.last.next = l2.first
    val first = l1.first
    val last = l2.last
    return new CoalescingList<T>(first, last)
  }
    
  def static <T> CoalescingList<T> build(T ... items) {
    var Node<T> first = null
    var Node<T> last = null
    var Node<T> prev = null
    for (var i = 0; i < items.size; i++) {
      val node = new Node<T>
      node.contents = items.get(i)
      if (prev !== null) {
        prev.next = node
      }
      if (i === 0) {
        first = node
      }
      if (i === items.size - 1) {
        last = node
      }
    }
    return new CoalescingList(first, last)
  }
  
  @Data
  static class CoalescingList<T> implements Iterable<T> {
    val Node<T> first
    val Node<T> last
    
    override iterator() {
      return new Iterator<T>() {
        var Node<T> current = new Node<T>() => [contents = null; next = first]
        override hasNext() {
          current !== last
        }
        override next() {
          current = current.next
          return current.contents
        }
      }
    }
  }
  
  static class Node<T> {
    T contents
    Node<T> next
  }
  
  @Data
  static class Island { 
    val String cell
    val String chr
    val Range<Integer> pos
    val Counter<Integer> posterior
    override toString() { '''Isl[«cell»,«chr»,«pos»,«posterior.totalCount»,«posterior.entries.entrySet.map[""+key+"="+value].join(",")»]''' }
  }
  
  def execDirs() { execDirs(directory, execInDirectory) }
  def static execDirs(File directory, String execInDirectory) {
    val result = BriefFiles.ls(directory).
      filter[isDirectory].
      map[new File(it, execInDirectory)].
      map[new ParsedExperiment(it)].
      toList
    Collections::sort(result, [x,y|cellId(x).compareTo(cellId(y))])
    return result
  }
  
  def static String cellId(ParsedExperiment parsed) {
    val cellPath = new File(parsed.arguments.get("model.data.readCounts.dataSource"))
    return cellPath.name.replaceAll("[.]csv", "")
  }
  
  static class Website extends BootstrapHTMLRenderer {
    val extension CombineCellPosteriors cellPost
    override htmlSupportFilesPrefix() { "../../../.html_support" }
    override container() { "container-fluid" }
    new(CombineCellPosteriors cellPost) {
      super(cellPost.directory.name)
      this.cellPost = cellPost
      documents.add(gather("Paths", [new File(it, "paths.pdf")]))
      documents.add(gather("Rates", [new File(it, "posteriorPlots/switchRate.pdf")]))
      documents.add(gather("Lambda", [new File(it, "monitoringPlots/lambdaInstantaneous.pdf")]))
    }
    
    def Document gather(String title, Function1<File,File> fileGrabber) {
      new Document(title) [
        for (execDir : execDirs) {
          val cellId = cellId(execDir)
          it += new Embed(fileGrabber.apply(execDir.execDir).absoluteFile) {
            override height() { "150px" }
            override title() { LINK(execDir.execDir.absolutePath) + "" + cellId + "(" + expectedNumberOfEvents.getCount(cellId) + " events)" + ENDLINK }
          }
        }
      ]
    }
  }
  
  def void callR(File _scriptFile, String commands) {
    val scriptFile = if (_scriptFile === null) BriefFiles.createTempFile() else _scriptFile
    BriefIO.write(scriptFile, commands);
    Command.call(Rscript.appendArg(scriptFile.getAbsolutePath()))
  }
  
  var Command _Rscript
  def Command Rscript() {
    if (_Rscript === null)
      _Rscript = Command.cmd(rCmd)
    return _Rscript
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
  
}