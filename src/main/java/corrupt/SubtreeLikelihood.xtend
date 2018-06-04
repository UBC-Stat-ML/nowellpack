package corrupt

import org.eclipse.xtend.lib.annotations.Accessors

/**
 * Product over all inclusion and exclusion probabilities of a subtree,
 * stored in log scale. 
 */
class SubtreeLikelihood {
  @Accessors(PUBLIC_GETTER) val double inclusionLog  // p
  @Accessors(PUBLIC_GETTER) val double exclusionLog  // q
  @Accessors(PUBLIC_GETTER) val double logProductPQ // Used for sampling: log ( \product_{children c} (p_c + q_c) 
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