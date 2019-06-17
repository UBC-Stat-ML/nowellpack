package humi.posets

import blang.inits.experiments.Experiment
import blang.inits.parsing.Posix
import humi.HumiStaticUtils
import java.io.File
import java.util.LinkedHashSet
import briefj.BriefIO
import humi.v5.DeltaMethod
import blang.inits.Arg
import blang.inits.DefaultValue
import bayonet.graphs.DotExporter

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
    plot(poset)
  }
  
  def File plot(Poset<LabelledInterval> poset) {
    val graph = GraphPoset.from(poset).graph
    val hasse = Posets.hasseDiagram(graph)
    val result = results.getFileInResultFolder("hasse.dot")
    new DotExporter(hasse).export(result) 
    return result
  }
  
  def loadIntervals() {
    val result = new LinkedHashSet<LabelledInterval>
    val names = new LinkedHashSet
    for (line : BriefIO::readLines(intervalsCSVFile).indexCSV) {
      val name = line.get(geneAField) + " (" + line.get(sgRNAField) + ")"
      if (names.contains(name)) throw new RuntimeException
      names.add(name)
      val left = Double.parseDouble(line.get(DeltaMethod.Columns::logRatioLeftBound.toString))
      val right = Double.parseDouble(line.get(DeltaMethod.Columns::logRatioRightBound.toString))
      val interval = new LabelledInterval(name, left, right)
      result.add(interval)
    }
    return result
  }
  
  def Poset<LabelledInterval> loadPoset() {
    if (threshold < 0.0) throw new RuntimeException
    val _objects = loadIntervals.filter[
      right < -threshold || left > threshold
    ].toSet
    _objects.add(new LabelledInterval("controls", 0.0, 0.0))
    return new Poset<LabelledInterval>() {
      override compare(LabelledInterval first, LabelledInterval second) {
        if (first.right < second.left) return -1
        if (second.right < first.left) return 1
        return null
      }
      override objects() {
        return _objects
      }
    }
  }
  
  def static void main(String [] args) {
    Experiment::start(args, Posix::parse(args), HumiStaticUtils::parsingConfigs)
  }
}