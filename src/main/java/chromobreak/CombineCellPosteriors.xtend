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
import java.util.List
import briefj.BriefMaps
import briefj.collections.Counter
import java.util.Collections
import com.google.common.collect.Comparators
import org.eclipse.xtext.xbase.lib.Functions.Function1

class CombineCellPosteriors extends Experiment {
  
  @Arg File directory
  
  @Arg @DefaultValue("results/latest")
  String execInDirectory = "results/latest"
  
  @Data
  static class Key {
    val String cell
    val String chr
    val int pos
  }
  
  override run() {
    new Website(this).renderInto(results.getFileInResultFolder("output"))
    // process posterior breakpoints
    // cell,chromo,pos,delta
    var changePoints = new LinkedHashMap<Key,Counter<Integer>>
    for (exec : execDirs) {
      val cellId = cellId(exec)
      val paths = new File(exec.execDir, "chromoplots/hmms.csv")
      var lines = BriefIO::readLines(paths).indexCSV.toList
      val nIters = lines.map[Integer.parseInt(get("sample"))].max
      for (i : 1 ..< lines.size - 1) {
        val state = lines.get(i).get("value")
        val pc = lines.get(i-1).get("map_key_0")
        val nc = lines.get(i+1).get("map_key_0")
        val sample = Integer.parseInt(lines.get(i).get("sample"))
        if (state == "NA" && pc == nc && sample > nIters/2) {
          val prev = Integer::parseInt(lines.get(i-1).get("value"))
          val next = Integer::parseInt(lines.get(i+1).get("value"))
          val delta = next - prev
          val key = new Key(cellId, nc, Integer.parseInt(lines.get(i).get("positions")))
          BriefMaps.getOrPut(changePoints, key, new Counter).incrementCount(delta, 1.0)
        }
      }
    }
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
  
  def execDirs() {
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
      documents.add(gather("Lambda", [new File(it, "monitoringPlots/lambdaInstantaneous.pdf")]))
    }
    
    def Document gather(String title, Function1<File,File> fileGrabber) {
      new Document(title) [
        for (execDir : execDirs) {
          val cellId = cellId(execDir)
          it += new Embed(fileGrabber.apply(execDir.execDir).absoluteFile) {
            override height() { "150px" }
            override title() { LINK(execDir.execDir.absolutePath) + "cell " + cellId + ENDLINK }
          }
        }
      ]
    }
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
  
}