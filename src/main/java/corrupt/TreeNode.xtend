package corrupt

import org.eclipse.xtend.lib.annotations.Data

/**
 * Either a locus, a special root node, or a cell. 
 */
@Data class TreeNode {
  val String id
  
  override toString() { id }
}
