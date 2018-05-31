package breaktree

import blang.mcmc.Samplers
import org.eclipse.xtend.lib.annotations.Data
import java.util.List
import java.util.ArrayList

@Samplers(SplitSampler)
@Data class Split {
  val Tree tree
  val int locus
  val List<TipIndicator> tipIndicators
  
  new (Tree tree, int locus, int nCells) {
    this.tree = tree
    this.locus = locus
    this.tipIndicators = new ArrayList(nCells)
    for (i : 0 ..< nCells)
      this.tipIndicators.add(new TipIndicator)
  }
  
  def nCells() {
    tipIndicators.size
  }
}
