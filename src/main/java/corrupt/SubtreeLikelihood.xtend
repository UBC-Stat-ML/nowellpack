package corrupt

import org.eclipse.xtend.lib.annotations.Data

/**
 * Product over all inclusion and exclusion probabilities of a subtree,
 * stored in log scale. 
 */
@Data class SubtreeLikelihood {
  val double inclusionLog  // p; in the notes: p_v^1
  val double exclusionLog  // q; in the notes: p_v^0
  val double logProductPQ  // log ( \product_{children c} (p_c + q_c) 
  new (double logP, double logQ, double logProductPQ) {
    inclusionLog = logP
    exclusionLog = logQ
    this.logProductPQ = logProductPQ
  }
  def static SubtreeLikelihood tip(double logP, double logQ) {
    return new SubtreeLikelihood(logP, logQ, 0.0)
  }
  def static SubtreeLikelihood missingTip() {
    return tip(0.0, 0.0)
  }
}