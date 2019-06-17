package humi.posets

import blang.inits.experiments.Experiment
import blang.inits.parsing.Posix
import humi.HumiStaticUtils
import blang.inits.Arg
import java.util.LinkedHashSet
import java.util.ArrayList
import briefj.collections.Counter
import java.util.Set
import java.util.LinkedHashMap
import static extension briefj.BriefMaps.getOrPutSet

class ComparePosets extends Experiment {
  
  @Arg 
  public Intervals2Poset condition1
  
  @Arg 
  public Intervals2Poset condition2
  
  override run() {
    condition1.results = results.child("condition1")
    condition2.results = results.child("condition2")
    
    condition1.run
    condition2.run
    
    val poset1 = condition1.loadPoset
    val poset2 = condition2.loadPoset
    
    val vertices = new LinkedHashSet(poset1.objects)
    vertices.retainAll(poset2.objects)
//    val vertexList = new ArrayList(vertices)
    
    val counter = new Counter<Pair<String,Boolean>> // 0: decrease in rank from c1 to c2; 1: increase
    val partners = new LinkedHashMap<Pair<String,Boolean>,Set<String>>
    
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
      
      println("" + v.key + " " + (if (v.value) "increased" else "decreased") + " in rank from condition 1 to 2 compared to:")
      
      for (p : ps) {
        println("\t" + p)
        val key = Pair.of(p, !v.value)
        counter.incrementCount(key, -1.0)
        partners.get(key).remove(v.key)
      }
    }
        
    
//    for (i : 0 ..< vertexList.size)
//       for (j : (i+1) ..< vertexList.size) {
//         val v1 = vertexList.get(i)
//         val v2 = vertexList.get(j)
//         if (reversed(poset1, poset2, v1, v2)) {
//           counter.incrementCount(v1, 1.0)
//           counter.incrementCount(v2, 1.0)
//         }
//       }
    

  }
  
  def static <T> boolean reversed(Poset<T> p1, Poset<T> p2, T v1, T v2) {
    val c1 = p1.compare(v1, v2)
    val c2 = p2.compare(v1, v2)
    if (c1 === null || c2 === null) return false
    if (c1 > 0 && c2 < 0) return true
    if (c2 > 0 && c1 < 0) return true
    return false
  }
  
  def static void main(String [] args) {
    Experiment::start(args, Posix::parse(args), HumiStaticUtils::parsingConfigs)
  }
}