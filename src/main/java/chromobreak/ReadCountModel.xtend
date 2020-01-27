package chromobreak

import org.eclipse.xtend.lib.annotations.Data
import blang.core.RealVar

@Data
class ReadCountModel {
  
  val RealVar a
  val RealVar h 
  val RealVar k 
  val RealVar sd 
  val RealVar sdSlope 
  
  def double logDensity(double logGC, double logReadCount, int state) {
    
    if (sd.doubleValue <= 0.0) blang.types.StaticUtils::invalidParameter
    
    if (Double.isNaN(logGC) || Double.isNaN(logReadCount))
      return 0.0 // skip missing observations 
    
    if (logReadCount == Double.NEGATIVE_INFINITY) {
      if (state === 0) return 0.0 else return Double.NEGATIVE_INFINITY
    }
    
    val lowerCutoff = mean(logGC, 0) - 3.0 * sd(logGC, 0)
    if (logReadCount < lowerCutoff && state === 0) {
      return Math.log(1.0/lowerCutoff)
    }
    val mean = mean(logGC, state)
    val sd = sd(logGC, state)
    if (sd <= 0.0) blang.types.StaticUtils::invalidParameter
    val result = (-0.5 * Math::pow( (logReadCount - mean) / sd, 2)) - Math.log(sd) -  LOG_SQRT_2_PI

    return result  
  }
  
  val static LOG_SQRT_2_PI = Math::log(Math::sqrt(2.0 * Math::PI))
  
  def double mean(double logGC, int state) {
    a.doubleValue * Math::pow(logGC - h.doubleValue, 2.0) + k.doubleValue + Math::log(state)
  }
  
  def double sd(double logGC, int state) { 
    sd.doubleValue + sdSlope.doubleValue * state
  }
}