package corrupt.distances

import corrupt.Locus
import corrupt.PerfectPhylo
import java.util.Set
import java.util.HashSet
import corrupt.TreeNode
import corrupt.DirectedTree

class LocusEdgeStat implements TreeStatistic<Pair<Locus,Locus>> {
  override compute(PerfectPhylo phylo) {
    locusEdges(phylo.collapsedTree)
  }
  
  static def Set<Pair<Locus,Locus>> locusEdges(DirectedTree<Set<TreeNode>> tree) {
    val result = new HashSet
    for (topNodeSet : tree.nodes) 
      for (bottomNodeSet : tree.children(topNodeSet))
        for (topNode : topNodeSet)
          if (topNode instanceof Locus)
            for (botNode : bottomNodeSet)
              if (botNode instanceof Locus)
                result.add((topNode as Locus) -> (botNode as Locus))
    return result
  }
}