package corrupt.distances

import corrupt.PerfectPhylo
import java.util.Set
import blang.inits.Implementations

@Implementations(LocusEdgeStat)
interface TreeStatistic<T> {
  def Set<T> compute(PerfectPhylo phylo)
}