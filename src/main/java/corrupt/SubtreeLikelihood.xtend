package corrupt

import org.eclipse.xtend.lib.annotations.Accessors

/**
 * Product over all inclusion and exclusion probabilities of a subtree,
 * stored in log scale. 
 */
class SubtreeLikelihood {
  @Accessors(PUBLIC_GETTER) val double inclusionLog 
  @Accessors(PUBLIC_GETTER) val double exclusionLog  
  new (double logP, double logQ) {
    inclusionLog = logP
    exclusionLog = logQ
  }
}