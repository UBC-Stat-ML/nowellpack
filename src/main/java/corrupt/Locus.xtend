package corrupt

import blang.inits.DesignatedConstructor
import blang.inits.Input

class Locus extends TreeNode {
  public static val PREFIX = "locus_"
  @DesignatedConstructor new(@Input String id) { super(id.replaceFirst(PREFIX,  "")) }
  override toString() { PREFIX + id }
}
