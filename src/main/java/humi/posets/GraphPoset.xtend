package humi.posets

import org.jgrapht.DirectedGraph
import org.eclipse.xtend.lib.annotations.Data
import bayonet.graphs.GraphUtils
import org.apache.commons.lang3.tuple.Pair

@Data
class GraphPoset<T> implements Poset<T> {
  
  val DirectedGraph<T, Pair<T,T>> graph
  
  /**
   * Return input if already a GraphPoset
   */
  def static <T> GraphPoset<T> from(Poset<T> poset) {
    if (poset instanceof GraphPoset) return poset as GraphPoset<T>
    val result = GraphUtils.newDirectedGraph
    for (v : poset.objects) result.addVertex(v)
    for (v1 : poset.objects)
      for (v2 : poset.objects)
        if (poset.compare(v1, v2) !== null && poset.compare(v1, v2) > 0)
          result.addEdge(v1, v2)
    return new GraphPoset(result)
  }
  
  override compare(T first, T second) {
    if (graph.containsEdge(first, second)) return 1
    if (graph.containsEdge(second, first)) return -1
    return null
  }
  override objects() {
    graph.vertexSet
  }
}