package corrupt

import blang.mcmc.Samplers
import org.eclipse.xtend.lib.annotations.Data
import java.util.List
import java.util.ArrayList
import static corrupt.Root.root

@Samplers(SplitSampler)
@Data class Split {
  val DirectedTree<TreeNode> tree
  val Locus locus
  val List<TipIndicator> tipIndicators
  
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
