package corrupt

import java.util.Map
import java.util.LinkedHashMap
import java.util.Set

import static corrupt.TreeNode.root
import org.eclipse.xtend.lib.annotations.Data

@Data class PerfectPhylo {
  val DirectedTree<TreeNode> tree
  val Map<Locus, Split> splits
  val Set<Locus> loci
  val Set<Cell> cells
  
  /**
   * Initialized with a star tree.
   */
  new(Set<Cell> cells, Set<Locus> loci) { 
    this.cells = cells
    this.loci = loci
    tree = new DirectedTree(root)
    splits = new LinkedHashMap
    for (cell : cells)
      tree.addEdge(root, cell)
    for (locus : loci) 
      splits.put(locus, Split::initializeEmpty(tree, locus, cells))
  }
  
  def TipIndicator tipIndicator(Cell cell, Locus locus) {
    return splits.get(locus).tipIndicators.get(cell)
  }
}