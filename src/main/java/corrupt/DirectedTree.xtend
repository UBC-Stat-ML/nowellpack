package corrupt

import java.util.List
import java.util.ArrayList
import org.eclipse.xtend.lib.annotations.Data
import org.eclipse.xtend.lib.annotations.Accessors
import static java.util.Collections.emptyList
import java.util.Collection
import bayonet.graphs.GraphUtils
import org.jgrapht.DirectedGraph

@Data class DirectedTree<T> { 
  @Accessors(NONE)
  val DirectedGraph<T, org.apache.commons.lang3.tuple.Pair<T,T>> graph
  val T root
  
  new(T root) {
    this(root,  GraphUtils.<T>newDirectedGraph)
    graph.addVertex(root)
  }
  
  new(T root, DirectedGraph<T, org.apache.commons.lang3.tuple.Pair<T,T>> graph) {
    this.graph = graph  
    this.root = root
  }
  
  def T parent(T node) {
    val iterator = graph.incomingEdgesOf(node).iterator
    if (!iterator.hasNext) return null
    val result = iterator.next
    if (iterator.hasNext) throw new RuntimeException
    return result.left
  }
  
  def List<T> children(T node) {
    val result = new ArrayList
    for (edge : graph.outgoingEdgesOf(node))
     result.add(edge.right)
    return result
  }
  
  def boolean isLeaf(T node) {
    return graph.outDegreeOf(node) === 0
  }
  
  def boolean isRoot(T node) {
    return graph.inDegreeOf(node) === 0
  }  
  
  def boolean hasNode(T node) {
    return graph.vertexSet.contains(node)
  }
  
  def Collection<T> nodes() { return graph.vertexSet }
  
  def List<T> collapseEdge(T bottomOfEdge) {
    val topOfEdge = parent(bottomOfEdge)
    if (topOfEdge === null) 
      throw new RuntimeException("This does not define an edge: no parent for " + bottomOfEdge)
    val movedChildren = children(bottomOfEdge)
    for (child : movedChildren) {
      graph.removeEdge(bottomOfEdge, child)
      graph.addEdge(topOfEdge, child)
    }
    graph.removeVertex(bottomOfEdge) 
    return movedChildren
  }
  
  def void addEdge(T existingTopNode, T newBottomNode) {
    addEdge(existingTopNode, newBottomNode, emptyList)
  }
  
  /**
   * See Figure in paper. The movedChildren are those that become the children 
   * of the newBottomNode (in bold in figure).  
   */
  def void addEdge(T existingTopNode, T newBottomNode, List<T> movedChildren) {
    if (!graph.addVertex(newBottomNode))
      throw new RuntimeException("Should not be in the tree: " + newBottomNode)
    if (!graph.vertexSet.contains(newBottomNode))
      throw new RuntimeException("Should be in the tree: " + existingTopNode) 
    graph.addEdge(existingTopNode, newBottomNode)
    for (child : movedChildren) {
      if (graph.removeEdge(existingTopNode, child) === null)
        throw new RuntimeException
      graph.addEdge(newBottomNode, child)
    }
  }
  
  override String toString() {
    return toString(0, root).toString
  }
  
  def private StringBuilder toString(int pad, T node) {
    val result = new StringBuilder
    result.append((0..<pad).map["  "].join(""))
    result.append(node + "\n")
    for (child : children(node)) 
      result.append(toString(pad+1, child))
    return result
  }
}
