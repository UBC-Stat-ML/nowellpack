package corrupt

import org.eclipse.xtend.lib.annotations.Data
import blang.inits.DesignatedConstructor
import blang.inits.Input

@Data class Cell extends TreeNode {
  public static val PREFIX = "cell_"
  override toString() { PREFIX + id }
  @DesignatedConstructor new(@Input String id) { 
  	super(id.replaceFirst(PREFIX,  ""))
  }
}