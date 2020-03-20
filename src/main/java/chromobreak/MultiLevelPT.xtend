package chromobreak

import blang.engines.internals.factories.PT
import blang.engines.internals.Spline.MonotoneCubicSpline
import java.util.List
import chromobreak.SingleCellHMMs.MultiLevel
import org.apache.commons.math3.analysis.UnivariateFunction
import blang.engines.internals.EngineStaticUtils
import java.util.ArrayList
import blang.inits.Arg
import blang.inits.DefaultValue

class MultiLevelPT extends PT {
  
  public static MultiLevel configs
  
  @Arg @DefaultValue("0.3")
  double minFractionPerLevel = 0.3
  
  override List<Double> fixedSizeOptimalPartition(MonotoneCubicSpline cumulativeLambdaEstimate, int size) {
    if (configs.n !== 2) throw new RuntimeException
    val lambda1 = cumulativeLambdaEstimate.interpolate(0.5)
    val lambda2 = cumulativeLambdaEstimate.interpolate(1.0) - lambda1
    val minNumberPerLevel = (size * minFractionPerLevel) as int
    
    // find optimal allocation 
    val n1 =findN1OutOf2(lambda1, lambda2, size, configs.b, minNumberPerLevel) 
    val n2 = size - n1 - 1
    
    val result = new ArrayList<Double>
    
    // sub-allocation within each computational regime
    var UnivariateFunction f = [x | cumulativeLambdaEstimate.interpolate(x / 2.0)]
    var list = EngineStaticUtils::fixedSizeOptimalPartition(f, n1 + 1)  // +1 b/c of avoid duplicate at boundary
    for (var i = 0 ; i < list.size - 1; i++) 
      result.add(list.get(i) / 2.0)
      
    f = [x | cumulativeLambdaEstimate.interpolate(0.5 + x / 2.0) - lambda1]
    list = EngineStaticUtils::fixedSizeOptimalPartition(f, n2 + 1) 
    for (var i = 0 ; i < list.size - 1; i++)
      result.add(0.5 + list.get(i) / 2.0)
      
    result.add(1.0)
    return result
  }
  
  static def int findN1OutOf2(double lambda1, double lambda2, int n, double discount, int minNumberPerLevel) {
    var min = Double.POSITIVE_INFINITY
    var argmin = -1
    for (k : minNumberPerLevel .. (n-minNumberPerLevel)) {
      val double n1 = k
      val double n2 = n - k
      val cur = 
        (
          lambda1/(1.0 - lambda1/n1) +
          lambda2/(1.0 - lambda2/n2) 
        ) * (n1 / discount + n2) 
      if (cur < min) {
        min = cur
        argmin = k
      }
    }
    return argmin
  }
}