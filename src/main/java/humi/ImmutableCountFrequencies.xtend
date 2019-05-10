package humi

import com.rits.cloning.Immutable
import java.util.Map
import java.util.LinkedHashMap

@Immutable
class ImmutableCountFrequencies implements CountFrequencies {
  val Map<Integer, Integer> data = new LinkedHashMap // count -> frequency
  val int nDataPoints
  
  new(Map<Integer, Integer> data) {
    this.data.putAll(data)
    this.nDataPoints = data.values.stream.mapToInt[intValue].sum
  }
  
  override distinctCounts() {
    return data.keySet
  }
  override frequency(int count) {
    return data.get(count).intValue
  }
  override nDataPoints() {
    return nDataPoints
  }
  
  override toString() { toString(this) }
}