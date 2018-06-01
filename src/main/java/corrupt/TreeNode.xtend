package corrupt

import org.eclipse.xtend.lib.annotations.Data

@Data class TreeNode {
  val int index
  val int type // using this to make sure hashcodes are nice
}
