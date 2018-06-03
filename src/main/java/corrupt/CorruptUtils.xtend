package corrupt

import java.util.Collection
import blang.types.Plate
import java.util.List
import bayonet.distributions.Random
import java.util.ArrayList

class CorruptUtils {
  def static Collection<Locus> set(Plate<String> lociDescriptions) {
    lociDescriptions.indices.map[new Locus(key)].toList 
  }
  def static Collection<Locus> syntheticLoci(int n) {
    (0 ..< n).map[new Locus("syntheticLocus_" + it)].toList
  }
  def static List<Locus> lociAndRoot(DirectedTree<TreeNode> tree) {
    tree.nodes.filter[it instanceof Locus].map[it as Locus].toList
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
