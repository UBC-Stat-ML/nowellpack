package humi

import java.util.Collection
import blang.inits.Input
import blang.inits.DesignatedConstructor
import java.util.LinkedHashMap

/**
 * Summarize list of counts by grouping identical counts together.
 * The frequency of a count is then the number of times it occurs in the list.
 * The number of data points is the number of items in the original list.
 */
interface CountFrequencies {
  def Collection<Integer> distinctCounts() 
  def int frequency(int count)
  def int nDataPoints()
  
  @DesignatedConstructor
  def static ImmutableCountFrequencies parse(
    @Input(formatDescription = "semi-column separated pairs each of format <count>x<frequency>, e.g. '3x4;2x1' means '4, 4, 4, 1, 1'") 
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
  
  def static String toString(CountFrequencies c) {
    return c.distinctCounts.map["" + c.frequency(it) + "x" + it].join(";")
  }
}