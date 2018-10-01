package humi

import com.rits.cloning.Immutable
import java.util.Map
import java.util.LinkedHashMap

@Immutable
class ImmutableCountFrequencies implements CountFrequencies {
  val Map<Integer, Integer> data = new LinkedHashMap
  
  new(Map<Integer, Integer> data) {
    data.putAll(data)
  }
  
  override counts() {
    return data.keySet
  }
  override frequency(int count) {
    return data.get(count).intValue
  }
  
}