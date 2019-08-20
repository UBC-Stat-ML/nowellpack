package humi

import blang.core.RealVar
import org.eclipse.xtend.lib.annotations.Data

@Data
class DerivedRealVar<T> implements RealVar {
  
  def static <T> DerivedRealVar<T> derivedReal(()=>Double function) {
    return new DerivedRealVar<T>(function)
  }
  
  val ()=>Double function
  override doubleValue() {
    return function.apply()
  }
  override toString() { return "" + doubleValue }
}