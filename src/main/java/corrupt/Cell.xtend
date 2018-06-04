package corrupt

import org.eclipse.xtend.lib.annotations.Data
import blang.inits.DesignatedConstructor
import blang.inits.Input

@Data class Cell extends TreeNode {
  override toString() { "cell_" + id }
  @DesignatedConstructor new(@Input String id) { super(id) }
}