package humi

import java.util.Collection
import blang.inits.Input
import blang.inits.DesignatedConstructor
import java.util.LinkedHashMap

interface CountFrequencies {
  def Collection<Integer> counts()
  def int frequency(int count)
  
  @DesignatedConstructor
  def static ImmutableCountFrequencies parse(
    @Input(formatDescription = "semi-column separated pairs each of format <count>x<frequency>, e.g. '3x0;2x1' means {3, 3, 3, 1, 1}") 
    String string
  ) {
    val data = new LinkedHashMap<Integer,Integer> => [
      for (pair : string.split("[;]")) {
        val splitPair = pair.split("[x]")
        if (splitPair.size !== 2) throw new RuntimeException
        put(Integer.parseInt(splitPair.get(0)), Integer.parseInt(splitPair.get(1)))
      }
    ]
    return new ImmutableCountFrequencies(data)
  }
}