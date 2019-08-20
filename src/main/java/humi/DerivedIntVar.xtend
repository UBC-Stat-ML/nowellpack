package humi

import org.eclipse.xtend.lib.annotations.Data
import blang.core.IntVar

@Data
class DerivedIntVar<T> implements IntVar {
  
  def static <T> DerivedIntVar<T> derivedInt(()=>Integer function) {
    return new DerivedIntVar<T>(function)
  }
  
  val ()=>Integer function
  override intValue() {
    return function.apply()
  }
  override toString() { return "" + intValue }
}