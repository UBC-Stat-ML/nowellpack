package corrupt

import blang.types.Plate
import java.util.List
import bayonet.distributions.Random
import java.util.ArrayList
import static corrupt.TreeNode.root
import java.util.Set

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
}
