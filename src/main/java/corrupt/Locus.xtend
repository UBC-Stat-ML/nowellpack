package corrupt

import org.eclipse.xtend.lib.annotations.Data

@Data class Locus implements TreeNode {
  public val static Locus root = new Locus("ROOT")
  val String id
  override toString() { id }
}
