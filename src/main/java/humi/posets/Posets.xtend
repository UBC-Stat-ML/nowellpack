package humi.posets

import org.jgrapht.DirectedGraph
import org.apache.commons.lang3.tuple.Pair
import bayonet.graphs.GraphUtils
import java.util.LinkedHashSet
import org.jgrapht.traverse.DepthFirstIterator
import java.util.Set
import com.google.common.collect.Sets
import bayonet.graphs.DotExporter
import java.io.File
import briefj.BriefIO
import humi.v5.DeltaMethod
import org.eclipse.xtend.lib.annotations.Data
import java.util.List
import java.util.Map
import briefj.BriefMaps
import java.util.LinkedHashMap

class Posets {
  
  def static <T> DirectedGraph<T, Pair<T,T>> toGraph(Poset<T> poset) {
    val result = GraphUtils.newDirectedGraph
    for (v : poset.objects) result.addVertex(v)
    for (v1 : poset.objects)
      for (v2 : poset.objects)
        if (poset.compare(v1, v2) !== null && poset.compare(v1, v2) > 0)
          result.addEdge(v1, v2)
    return result
  }
  
  def static <T> DirectedGraph<T, Pair<T,T>> hasseDiagram(DirectedGraph<T, Pair<T,T>> poset) {
    val result = GraphUtils.newDirectedGraph()
    for (v : poset.vertexSet) result.addVertex(v)
    for (v : poset.vertexSet) {
      val reducedChildren = new LinkedHashSet(poset.outgoingEdgesOf(v).map[value].toSet)
      for (child : new LinkedHashSet(reducedChildren)) {
        for (visited : [new DepthFirstIterator(poset, child)]) {
          if (visited != child && reducedChildren.contains(visited))
            reducedChildren.remove(visited) 
        }
      }
      for (c : reducedChildren)
        result.addEdge(v, c)
    }
    return result
  }
  
  def static Poset<LabelledInterval> fromIntervals(File csv) {
    val _objects = new LinkedHashSet<LabelledInterval>
    for (line : BriefIO.readLines(csv).indexCSV) {
      val name = line.get("gene") + " (" + line.get("sgrna") + ")"
      val left = Double.parseDouble(line.get(DeltaMethod.Columns::logRatioLeftBound.toString))
      val right = Double.parseDouble(line.get(DeltaMethod.Columns::logRatioRightBound.toString))
      val interval = new LabelledInterval(name, left, right)
      if (right < -0.5 || left > 0.5)
        _objects.add(interval)
    }
    _objects.add(new LabelledInterval("reference", 0.0, 0.0))
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
  
  def static Poset<Set<Integer>> power(int size) {
    val refSet = new LinkedHashSet<Integer> => [for (i : 0 ..< size) add(i)]
    val _objects = Sets.powerSet(refSet)
    return new Poset<Set<Integer>>() {
      override compare(Set<Integer> first, Set<Integer> second) {
        if (first == second) return 0
        if (first.containsAll(second)) return 1
        if (second.containsAll(first)) return -1
        return null
      }
      override objects() {
        return _objects
      }
    }
  }
  
  
  def static <T> DirectedGraph<List<T>, Pair<List<T>,List<T>>> simplify(DirectedGraph<T, Pair<T,T>> poset) {
    val eqClasses = new LinkedHashMap<Pair<Set<Pair<T,T>>, Set<Pair<T,T>>>, List<T>>
    val labels = new LinkedHashMap<T, List<T>>
    for (v : poset.vertexSet) {
      val key = Pair.of(poset.incomingEdgesOf(v), poset.outgoingEdgesOf(v))
      val list = BriefMaps::getOrPutList(eqClasses, key)
      list.add(v)
      labels.put(v, list)
    }
    val DirectedGraph<List<T>, Pair<List<T>,List<T>>> result = GraphUtils.newDirectedGraph()
    for (value : eqClasses.values) result.addVertex(value)
    for (value : eqClasses.values) {
      val rep = value.get(0)
      for (child : poset.outgoingEdgesOf(rep)) 
        result.addEdge(value, labels.get(child.right)) 
    }
    return result
  }
  
  def static void main(String [] args) {
    val poset = fromIntervals(new File("/Users/bouchard/experiments/humi/results/all/2019-06-16-22-04-03-vjVFcWrp.exec/estimates-subset2.csv"))
    val graph = /*simplify*/(toGraph(poset))
    new DotExporter(graph).export(new File("test.dot"))
    val hasse = hasseDiagram(graph)
    new DotExporter(hasse).export(new File("hasse.dot")) 
  }
}