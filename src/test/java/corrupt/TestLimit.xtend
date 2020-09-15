package corrupt

import blang.types.StaticUtils

class TestLimit {
  def static void main(String [] args) {
    val int n = 56
    val double density = 0.05
    
    val estimate = StaticUtils::logFactorial(n) + Math::log(1.0 - density) * (40*56) //n * n * (1.0 - Math::PI / 4.0)
    println(estimate)
    
    // n! * Pr(first col ok) * Pr(second col ok) ...
    
    var sum = StaticUtils::logFactorial(n)
    
    for (i : 1 .. n) {
      val int numberInCol = (density * n) as int
      sum += numberInCol * Math::log(i as double / n)
    }
    
    println(sum)
  }
}