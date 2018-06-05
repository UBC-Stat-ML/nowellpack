package corrupt

import java.util.Set
import blang.types.Plate
import java.util.List

import static corrupt.CorruptStaticUtils.*
import java.util.ArrayList
import bayonet.distributions.Random
import java.util.Collection

import static bayonet.graphs.GraphUtils.newUndirectedGraph
import java.util.LinkedHashSet

class CorruptExtensionUtils {
  
  def static <T> Set<T> set(Plate<T> plate) {
    plate.indices.map[key].toSet
  }
  
  def static List<TreeNode> lociAndRoot(DirectedTree<TreeNode> tree) {
    tree.nodes.filter[it instanceof Locus || it == root].toList
  }
  
  def static List<Cell> cells(DirectedTree<TreeNode> tree) {
    tree.nodes.filter[it instanceof Cell].map[it as Cell].toList
  }
  
  def static <T> T uniformElement(Random rand, List<T> list) {
    return list.get(rand.nextInt(list.size))
  }
   
  def static <T> List<T> uniformSubset(Random rand, List<T> list) {
    val result = new ArrayList
    for (item : list)
      if (rand.nextBernoulli(0.5))
        result.add(item)
    return result
  }
  
  def static uniformUndirectedTree(Random rand, int size) {
    val result = newUndirectedGraph
    if (size < 0) throw new RuntimeException
    if (size == 0) return result
    val free = new LinkedHashSet((0 ..< size).toSet)
    val pointers = newIntArrayOfSize(size)
    val firstNode = rand.nextInt(size)
    free.remove(firstNode) 
    pointers.set(firstNode, inTree)
    result.addVertex(firstNode)
    while (result.vertexSet.size < size) {
      val branch = loopErasedWalk(rand, pointers, rand.uniformElement(free))
      free.removeAll(branch)
      for (i : 1 ..< branch.size) {
        result.addVertex(branch.get(i-1))
        result.addVertex(branch.get(i))
        result.addEdge(branch.get(i-1), branch.get(i))
      }
    }
    return result
  }
  
  static private val int inTree = -1
  def private static List<Integer> loopErasedWalk(
    Random rand,
    int [] pointers, 
    int startPoint
  ) {
    val size = pointers.size
    
    // do the walk
    var current = startPoint
    while (pointers.get(current) !== inTree) {
      val next = rand.nextInt(size)
      pointers.set(current, next)
      current = next
    }
    
    // create the list
    var result = new ArrayList
    current = startPoint
    result.add(current)
    while (pointers.get(current) !== inTree) {
      current = pointers.get(current)
      result.add(current)
    }
    
    // reset all indices
    for (i : 0 ..< size)
      if (pointers.get(i) !== inTree)
        pointers.set(i, 0)
    
    // inscribe the tree
    for (entry : result)
      pointers.set(entry, inTree)
    
    return result
  }
  
  def private static <T> T uniformElement(Random rand, Collection<T> items) { 
    var index = rand.nextInt(items.size)
    for (item : items)
      if (index-- === 0)
        return item
    throw new RuntimeException
  }
}