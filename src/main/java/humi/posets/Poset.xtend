package humi.posets

import java.util.Collection

interface Poset<T> {
  
  /**
   * null if not comparable
   */
  def Integer compare(T first, T second)
 
 
  def Collection<T> objects()
}