package corrupt

import blang.types.Plate
import java.util.List
import bayonet.distributions.Random
import java.util.ArrayList
import static corrupt.TreeNode.root
import java.util.Set
import briefj.collections.UnorderedPair
import bayonet.math.SamplingUtils
import java.util.Collection
import java.util.Queue
import java.util.Collections
import org.jgrapht.Graph
import bayonet.graphs.GraphUtils

import static bayonet.graphs.GraphUtils.newUndirectedGraph
import org.jgrapht.UndirectedGraph
import java.util.LinkedHashSet

class CorruptUtils {
  
  def static <T> Set<T> set(Plate<T> plate) {
    plate.indices.map[key].toSet
  }
  
  def static Set<Locus> syntheticLoci(int n) {
    (0 ..< n).map[new Locus(it.toString)].toSet
  }
  
  def static Set<Cell> syntheticCells(int n) {
    (0 ..< n).map[new Cell(it.toString)].toSet
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
  
  def static double logNPerfectPhylo(int nCells, int nLoci) {
    return (nLoci + nCells - 1) * Math.log(nLoci + 1)
  }
  
  def static PerfectPhylo sampleUniformPerfectPhylo(Random rand, Set<Cell> cells, Set<Locus> loci) {
    val result = new PerfectPhylo(cells, loci)
    result.sampleUniform(rand)
    return result
  }
  
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
  
  static private val int inTree = -1
  def static sampleUniform(Random rand, int size) {
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
  
  def private static <T> T uniformElement(Random rand, Collection<T> items) {
    var index = rand.nextInt(items.size)
    for (item : items)
      if (index-- == 0)
        return item
    throw new RuntimeException
  }
  
  
//  def static UndirectedGraph<Integer,UnorderedPair<Integer,Integer>> sampleUniformUnrooted(Random rand, int nNodes)
//  {
//    val permutation = new ArrayList((0 ..< nNodes).toList)
//    Collections.shuffle(permutation, rand) 
//    val result = newUndirectedGraph
//    if (nNodes < 0) throw new RuntimeException
//    if (nNodes == 0) return result
//    result.addVertex(permutation.get(0))
//    while (result.vertexSet.size < nNodes) {
//      // start random walk
//      val treeSize = result.vertexSet.size
//      var prev = rand.nextInt(treeSize) 
//      var walkSize = 0
//      var boolean stop
//      do {
//        walkSize++
//        val next = result.vertexSet.size
//        result.addVertex(permutation.get(next))
//        result.addEdge(permutation.get(next), permutation.get(prev))
//        prev = next
//        stop = rand.nextBernoulli(treeSize as double / (nNodes as double - walkSize as double))
//      } while (!stop)
//    }
//    return result
//  }
}
