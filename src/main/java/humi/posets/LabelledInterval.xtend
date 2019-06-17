package humi.posets

import org.eclipse.xtend.lib.annotations.Data

@Data
class LabelledInterval {
  val String name
  val double left
  val double right
  override toString() { name }
}