package corrupt

import org.eclipse.xtend.lib.annotations.Data

/**
 * Either a locus, a special root node, or a cell. 
 */
@Data class TreeNode {
  val String id
  public static val TreeNode root = new TreeNode("ROOT")
  override toString() { id }
}
