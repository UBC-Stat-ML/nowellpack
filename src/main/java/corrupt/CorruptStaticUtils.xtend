package corrupt

import blang.types.Plate
import java.util.List
import bayonet.distributions.Random
import java.util.ArrayList
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

class CorruptStaticUtils {
  
  public static val TreeNode root = new TreeNode("ROOT")
  
  def static Set<Locus> syntheticLoci(int n) {
    (0 ..< n).map[new Locus(it.toString)].toSet
  }
  
  def static Set<Cell> syntheticCells(int n) {
    (0 ..< n).map[new Cell(it.toString)].toSet
  }
  
  def static double logNPerfectPhylo(int nCells, int nLoci) {
    return (nLoci + nCells - 1) * Math.log(nLoci + 1)
  }

}
