package humi.posets

import java.util.Collection

interface Poset<T> {
  
  /**
   * Same as java.util.Comparator but allows
   * null if the two objects are not comparable
   */
  def Integer compare(T first, T second)
 
 
  /**
   * All objects allowed to be compared.
   */
  def Collection<T> objects()
}