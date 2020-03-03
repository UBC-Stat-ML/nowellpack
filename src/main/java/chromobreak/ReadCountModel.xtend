package chromobreak

import org.eclipse.xtend.lib.annotations.Data
import blang.core.RealVar
import blang.inits.experiments.tabwriters.TidilySerializable
import blang.inits.experiments.tabwriters.TidySerializer.Context
import java.util.LinkedHashMap

@Data
class ReadCountModel implements TidilySerializable {
  
  public val static double x0 = -1.0
  val RealVar f0
  val RealVar f1
  val RealVar f2 
  val RealVar sd 
  val RealVar sdSlope 
  
  def double logDensity(double logGC, double logReadCount, int state) {
    
    if (sd.doubleValue <= 0.0) 
      blang.types.StaticUtils::invalidParameter
    
    if (Double.isNaN(logGC) || Double.isNaN(logReadCount))
      return 0.0 // treat as missing 
    
    if (state === 0) { // loss state
      if (logReadCount === Double.NEGATIVE_INFINITY)
        return 0.0
    
      val bound = mean(logGC, 1)
      if (!(bound > 0.0))
        return Double.NEGATIVE_INFINITY
      if (logReadCount > bound) 
        return Double.NEGATIVE_INFINITY
      else
        return -Math::log(bound)
    }
    
    if (logReadCount == Double.NEGATIVE_INFINITY) {
      return Double.NEGATIVE_INFINITY
    }
    
    val mean = mean(logGC, state)
    val sd = sd(logGC, state)
    if (sd <= 0.0) 
      blang.types.StaticUtils::invalidParameter
    return -0.5 * Math::pow( (logReadCount - mean) / sd, 2) - Math.log(sd) - CONSTANT
  }
  
  val static double CONSTANT = Math::log(Math::sqrt(2.0 * Math::PI))
  
  def double mean(double logGC, int state) {
    mean(logGC, state, f0.doubleValue, f1.doubleValue, f2.doubleValue)
  }
  
  def static double mean(double logGC, int state, double f0, double f1, double f2) {
    val a = f2 / 2.0
    val b = f1 - 2 * a * x0
    val c = f0 - a * x0 * x0 - b * x0
    a * logGC * logGC + b * logGC + c + Math::log(state)
  }
  
  def double sd(double logGC, int state) { 
    sd.doubleValue + sdSlope.doubleValue * state
  }
  
  override serialize(Context context) {
    for (state : 1 .. ChromoPostProcessor::nStates) {
      val evaluations = new LinkedHashMap<Double, Double>
      for (var double lgc = -1.1; lgc < -0.5; lgc += 0.05)
        evaluations.put(lgc, mean(lgc, state))
      context.recurse(evaluations, "state", state)
    }
  } 
}