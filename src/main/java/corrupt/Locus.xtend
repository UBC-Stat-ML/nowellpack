package corrupt

import blang.inits.DesignatedConstructor
import blang.inits.Input

class Locus extends TreeNode {
  @DesignatedConstructor new(@Input String id) { super(id) }
  override toString() { "locus_" + id }
}
