package corrupt

import java.util.Map
import java.util.LinkedHashMap
import java.util.Set

import static extension corrupt.CorruptExtensionUtils.*
import static corrupt.CorruptStaticUtils.*
import org.eclipse.xtend.lib.annotations.Data
import bayonet.distributions.Random
import briefj.BriefLog

@Data class PerfectPhylo {
  // Warning: updates on the tree need to be mirrored to the splits
  val DirectedTree<TreeNode> tree 
  val Map<Locus, Split> splits
  val Set<Cell> cells
  
  /**
   * Initialized with a star tree.
   */
  new(Set<Cell> cells, Set<Locus> loci) { 
    this.cells = cells
    tree = new DirectedTree(root)
    splits = new LinkedHashMap
    for (cell : cells)
      tree.addEdge(root, cell)
    for (locus : loci) 
      splits.put(locus, Split::initializeEmpty(tree, locus, cells))
  }
  
  def Set<Locus> getLoci() { splits.keySet }
  
  def TipIndicator tipIndicator(Cell cell, Locus locus) {
    return splits.get(locus).tipIndicators.get(cell)
  }
  
  def void updateAllSplits() {
    for (split : splits.values)
      split.updateTips
  }
  
  def void sampleUniform(Random rand) {
    BriefLog.warnOnce("Sample Uniform not  yet implemented!")
//    for (node : tree.nodes)
//      if (node !== root)
//        tree.collapseEdge(node)
//    val lociTopology = rand.sampleUniformUndirectedTree(loci.size + 1)
    updateAllSplits
  }
  
  override String toString() { tree.toString }
}