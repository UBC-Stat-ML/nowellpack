package corrupt

import blang.mcmc.Samplers
import java.util.List
import java.util.ArrayList
import static corrupt.Root.root
import java.util.Collection
import org.eclipse.xtend.lib.annotations.Accessors

@Samplers(SplitSampler)
class Split {
  val DirectedTree<TreeNode> tree
  val Locus locus
  @Accessors(PACKAGE_GETTER)
  val List<TipIndicator> tipIndicators
  
  def void moveTo(TreeNode parent, List<TreeNode> movedChildren) {
    tree.collapseEdge(locus)
    tree.addEdge(parent, locus, movedChildren)
    updateTips
  }
  
  def Collection<TreeNode> children(TreeNode node) {
    tree.children(node)
  }
  
  def private void updateTips() { _updateTips(tree.root, false) }
  def private void _updateTips(TreeNode node, boolean active) {
    switch node {
      Cell : {
        if (!tree.isLeaf(node)) throw new RuntimeException
        tip(node).included = active
      }
      default : {
        val childrenActive = active || node == locus
        for (child : tree.children(node)) 
          _updateTips(child, childrenActive)
      }
    }
  }
  
  def TipIndicator tip(Cell cell) {
    return tipIndicators.get(cell.index)
  }
  
  /**
   * Initialize the split to a leaf connected to the root. 
   * I.e. no cell having the trait. 
   */
  new (DirectedTree<TreeNode> tree, Locus locus, int nCells) {
    this.tree = tree
    if (tree.hasNode(locus))
      throw new RuntimeException
    this.locus = locus
    this.tipIndicators = new ArrayList(nCells)
    for (i : 0 ..< nCells)
      this.tipIndicators.add(new TipIndicator)
    tree.addEdge(root, locus)
  }
  
  def nCells() {
    tipIndicators.size
  }
}
