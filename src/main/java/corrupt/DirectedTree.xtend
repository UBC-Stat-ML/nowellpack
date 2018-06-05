package corrupt

import java.util.List
import java.util.ArrayList
import org.eclipse.xtend.lib.annotations.Data
import org.eclipse.xtend.lib.annotations.Accessors
import static java.util.Collections.emptyList
import java.util.Collection
import bayonet.graphs.GraphUtils
import java.util.Map
import java.util.Set
import java.util.LinkedHashMap
import java.util.LinkedHashSet

@Data class DirectedTree<T> { 
  val Map<T, Set<T>> childrenPtrs = new LinkedHashMap
  val Map<T, T> parentPtrs = new LinkedHashMap
  val T root
  
  new(T root) {
    this.root = root
    childrenPtrs.put(root, new LinkedHashSet)
  }
  
  def T parent(T node) {
    if (!hasNode(node))
      throw new RuntimeException
    return parentPtrs.get(node)
  }
  
  def Collection<T> children(T node) {
    if (!hasNode(node))
      throw new RuntimeException
    return childrenPtrs.get(node)
  }
  
  def boolean isLeaf(T node) {
    return children(node).size === 0
  }
  
  def boolean hasNode(T node) {
    return nodes.contains(node)
  }
  
  def Collection<T> nodes() { return childrenPtrs.keySet }
  
  def List<T> collapseEdge(T bottomOfEdge) {
    val topOfEdge = parent(bottomOfEdge)
    if (topOfEdge === null) 
      throw new RuntimeException("This does not define an edge: no parent for " + bottomOfEdge)
    val movedChildren = children(bottomOfEdge)
    for (child : movedChildren) {
      parentPtrs.put(child, topOfEdge)
      childrenPtrs.get(topOfEdge).add(child)
    }
    parentPtrs.remove(bottomOfEdge)
    childrenPtrs.remove(bottomOfEdge)
    childrenPtrs.get(topOfEdge).remove(bottomOfEdge)
    return new ArrayList(movedChildren)
  }
  
  def void addEdge(T existingTopNode, T newBottomNode) {
    addEdge(existingTopNode, newBottomNode, emptyList)
  }

  def void addEdge(T existingTopNode, T newBottomNode, List<T> movedChildren) {
    if (hasNode(newBottomNode))
      throw new RuntimeException("Should not be in the tree: " + newBottomNode)
    if (!hasNode(existingTopNode))
      throw new RuntimeException("Should be in the tree: " + existingTopNode) 
    childrenPtrs.put(newBottomNode, new LinkedHashSet)
    parentPtrs.put(newBottomNode, existingTopNode)
    childrenPtrs.get(existingTopNode).add(newBottomNode)
    for (child : movedChildren) {
      parentPtrs.put(child, newBottomNode)
      childrenPtrs.get(existingTopNode).remove(child)
      childrenPtrs.get(newBottomNode).add(child)
//      if (graph.removeEdge(existingTopNode, child) === null)
//        throw new RuntimeException
//      graph.addEdge(newBottomNode, child)
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
