package corrupt

import java.util.List
import static corrupt.Root.root
import java.util.ArrayList

class PerfectPhylo {
  val DirectedTree<TreeNode> tree
  val List<Split> splits
  
  // TODO: might be better to initialize with all mutations above? propose static methods
  
  /**
   * Initialized with a star tree.
   */
  new(int nCells, int nLoci) { 
    tree = new DirectedTree(root)
    splits = new ArrayList
    for (cellIndex : 0 ..< nCells)
      tree.addEdge(root, new Cell(cellIndex))
    for (locusIndex : 0 ..< nLoci) 
      splits.add(new Split(tree, new Locus(locusIndex), nCells))
  }
  
  def TipIndicator tipIndicator(int cellIndex, int locusIndex) {
    return splits.get(locusIndex).tipIndicators.get(cellIndex)
  }
}