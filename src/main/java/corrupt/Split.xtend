package corrupt

import blang.mcmc.Samplers
import java.util.LinkedHashMap
import java.util.Map
import java.util.Set
import org.eclipse.xtend.lib.annotations.Data

import static extension corrupt.CorruptExtensionUtils.*
import static corrupt.CorruptStaticUtils.*

@Samplers(PerfectPhyloGibbsSampler)
@Data class Split {
  // Warning: updates on the tree need to be mirrored to the splits
  val DirectedTree<TreeNode> tree
  val Locus locus
  val Map<Cell, TipIndicator> tipIndicators
  
  def public void updateTips() { _updateTips(tree.root, false) }
  def private void _updateTips(TreeNode node, boolean active) {
    switch node {
      Cell : {
        if (!tree.isLeaf(node)) throw new RuntimeException
        tipIndicators.get(node).included = active
      }
      default : {
        val childrenActive = active || node == locus
        for (child : tree.children(node)) 
          _updateTips(child, childrenActive)
      }
    }
  }
  
  /**
   * Initialize the split to a leaf connected to the root. 
   * I.e. no cell having the corresponding trait. 
   */
  def static Split initializeEmpty(DirectedTree<TreeNode> tree, Locus locus, Set<Cell> cells) {
    val Map<Cell, TipIndicator> tips = new LinkedHashMap
    for (cell : cells) 
      tips.put(cell, new TipIndicator(cell))
    return new Split(tree, locus, tips) 
  }
}
