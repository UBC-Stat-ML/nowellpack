package humi.posets

import org.jgrapht.DirectedGraph
import org.apache.commons.lang3.tuple.Pair
import bayonet.graphs.GraphUtils
import java.util.LinkedHashSet
import org.jgrapht.traverse.DepthFirstIterator
import java.util.Set
import com.google.common.collect.Sets
import java.util.List
import briefj.BriefMaps
import java.util.LinkedHashMap
import bayonet.graphs.DotExporter
import java.io.File

class Posets {
  
  def static dotExporter(Poset<String> poset) {
    val graph = GraphPoset.from(poset).graph
    val hasse = Posets.hasseDiagram(graph)
    return new DotExporter(hasse)
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
    val poset = power(3)
    val fullGraph = GraphPoset::from(poset).graph
    new DotExporter(fullGraph).export(new File("full.dot"))
    new DotExporter(hasseDiagram(fullGraph)).export(new File("hasse.dot"))
  }
}