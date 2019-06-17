package humi.posets

import blang.inits.experiments.Experiment
import blang.inits.parsing.Posix
import humi.HumiStaticUtils
import blang.inits.Arg
import java.util.LinkedHashSet
import briefj.collections.Counter
import java.util.Set
import java.util.LinkedHashMap
import static extension briefj.BriefMaps.getOrPutSet
import blang.inits.DefaultValue

class ComparePosets extends Experiment {
  
  @Arg 
  public Intervals2Poset condition1
  
  @Arg 
  public Intervals2Poset condition2
  
  @Arg              @DefaultValue("condition1")
  public String condition1Label = "condition1"
  
  @Arg              @DefaultValue("condition2")
  public String condition2Label = "condition2"
  
  override run() {
    val poset1 = condition1.loadPoset
    val poset2 = condition2.loadPoset
    
    comparePosets(poset1, poset2)
    
    for (entry : #[condition1Label -> poset1, condition2Label -> poset2]) {
      val exporter = Posets::dotExporter(entry.value)
      exporter.addVertexAttribute("fontcolor", [
        if (increases.containsKey(it) && decreases.containsKey(it)) return "yellow"
        else if (increases.containsKey(it)) return if (entry.key == condition2Label) "green" else "red"
        else if (decreases.containsKey(it)) return if (entry.key == condition2Label) "red" else "green"
        else return "black"
      ])
      exporter.export(results.getFileInResultFolder(entry.key + ".dot"))
    }
  }
  
  // non-redundant record of the increases and decreases
  public val increases = new LinkedHashMap
  public val decreases = new LinkedHashMap
  
  def <T> comparePosets(Poset<T> poset1, Poset<T> poset2) {
    val vertices = new LinkedHashSet(poset1.objects)
    vertices.retainAll(poset2.objects)
    
    val counter = new Counter<Pair<T,Boolean>> // 0: decrease in rank from c1 to c2; 1: increase
    val partners = new LinkedHashMap<Pair<T,Boolean>,Set<T>>
    
    for (v1 : vertices) 
      for (v2 : vertices) {
        val c1 = poset1.compare(v1, v2)
        val c2 = poset2.compare(v1, v2)
        if (c1 !== null && c2 !== null) {
          if (c1 == -1 && c2 == 1) { // symmetric case handled by loop structure
            /* record increase: */ { val key = Pair.of(v1, true); counter.incrementCount(key, 1.0); partners.getOrPutSet(key).add(v2) }            
            /* record decrease: */ { val key = Pair.of(v2, false);counter.incrementCount(key, 1.0); partners.getOrPutSet(key).add(v1) }
          }
        }
      }
    
    while (counter.totalCount > 0.0) {
      val v = counter.argMax
      counter.setCount(v, 0.0)
      
      val ps = partners.get(v)
      partners.put(v, null)
      
      println("" + v.key + " " + (if (v.value) "increased" else "decreased") + " in rank from " + condition1Label + " to " + condition2Label + " compared to:")
      
      (if (v.value) increases else decreases).put(v.key, new LinkedHashSet(ps))
      
      for (p : ps) {
        println("\t" + p)
        val key = Pair.of(p, !v.value)
        counter.incrementCount(key, -1.0)
        partners.get(key).remove(v.key)
      }
    }
  }
  
  def static void main(String [] args) {
    Experiment::start(args, Posix::parse(args), HumiStaticUtils::parsingConfigs)
  }
}