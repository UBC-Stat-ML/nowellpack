package corrupt

import static corrupt.Locus.root
import java.util.Collection
import java.util.Map
import java.util.LinkedHashMap

class PerfectPhylo {
  val protected DirectedTree<TreeNode> tree
  val protected Map<Locus, Split> splits
  val int nCells
  
  /**
   * Initialized with a star tree.
   */
  new(int nCells, Collection<Locus> loci) { 
    this.nCells = nCells
    tree = new DirectedTree(root)
    splits = new LinkedHashMap
    for (cellIndex : 0 ..< nCells)
      tree.addEdge(root, new Cell(cellIndex))
    for (locus : loci) 
      splits.put(locus, new Split(tree, locus, nCells))
  }
  
  def int nCells() { return nCells }
  def int nLoci() { return splits.size }
  
  def TipIndicator tipIndicator(Cell cell, Locus locus) {
    return splits.get(locus).tip(cell)
  }
  
  override toString() { tree.toString }
}