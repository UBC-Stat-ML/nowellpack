package humi.posets

import org.eclipse.xtend.lib.annotations.Data

@Data
class LabelledInterval {
  val String name   // should be unique
  val double left   // ignored in equal/hashcode (needed for ComparePosets)
  val double right  // ignored in equal/hashcode (needed for ComparePosets)
  override toString() { name }
  
  override hashCode() {
    name.hashCode
  }
  
  override equals(Object obj) {
    if (this === obj)
      return true;
    if (obj === null)
      return false;
    if (getClass() != obj.getClass())
      return false;
    val other = obj as LabelledInterval
    return name == other.name
  }
}