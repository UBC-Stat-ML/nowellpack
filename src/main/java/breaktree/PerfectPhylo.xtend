package breaktree

import org.eclipse.xtend.lib.annotations.Data

@Data class PerfectPhylo {
  val Tree tree
  val Split [] splits
  
  new(int nCells, int nLoci) { 
    tree = new Tree(nCells, nLoci)
    splits = newArrayOfSize(nLoci)
    for (locus : 0 ..< nLoci) 
      splits.set(locus, new Split(tree, locus, nCells))
  }
  
  def TipIndicator tipIndicator(int cellIndex, int locusIndex) {
    return splits.get(locusIndex).tipIndicators.get(cellIndex)
  }
}