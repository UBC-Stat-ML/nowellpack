package corrupt.post

import org.eclipse.xtend.lib.annotations.Data
import corrupt.Cell
import corrupt.Locus
import blang.core.RealVar
import blang.inits.ConstructorArg
import blang.inits.DesignatedConstructor
import blang.core.RealConstant

@Data class NoisyBinaryCLMatrix implements CellLocusMatrix {
  val ReadOnlyCLMatrix binaryMatrix 
  val RealVar fpRate
  val RealVar fnRate
  val int nPositive
  val int nNegative
  
  @DesignatedConstructor
  new(
    @ConstructorArg("binaryMatrix") ReadOnlyCLMatrix binaryMatrix, 
    @ConstructorArg("fpRate") double fpRate, 
    @ConstructorArg("fnRate") double fnRate
  ) {
    this(binaryMatrix, new RealConstant(fpRate), new RealConstant(fnRate))
  }
  
  new(
    ReadOnlyCLMatrix binaryMatrix, 
    RealVar fpRate, 
    RealVar fnRate
  ) {
    this.binaryMatrix = binaryMatrix
    this.fpRate = fpRate
    this.fnRate = fnRate 
    var int nPositive = 0
    var int nNegative = 0
    for (locus : binaryMatrix.loci)
      for (cell : binaryMatrix.cells)
        if (binaryValue(cell, locus))
          nPositive++
        else
          nNegative++
    this.nPositive = nPositive
    this.nNegative = nNegative
  }
  
  def boolean binaryValue(Cell cell, Locus locus) {
    switch(binaryMatrix.get(cell, locus)) {
        case 0.0 : false
        case 1.0 : true
        default  : 
          throw new RuntimeException
    }    
  }
  
  override get(Cell cell, Locus locus) { 
    val boolValue = binaryValue(cell, locus)
    check
    if (boolValue) return p_x1_given_y1
    else           return p_x1_given_y0
  }
  
  def boolean valid() {
    return fp > 0.0 && fp < 1.0 && fn > 0.0 && fn < 1.0
  }
  
  def void check() {
    if (!valid) 
      throw new RuntimeException
  }
  
  def double fp() { fpRate.doubleValue }
  def double fn() { fnRate.doubleValue }
  
  /**
   * P(X = 1 | Y = 1)
   */
  def double p_x1_given_y1() {
    return (1.0 - fn) / (1.0 - fn + fp)
  }
  
  /**
   * P(X = 1 | Y = 0)
   */
  def double p_x1_given_y0() {
    return fn / (fn + 1.0 - fp)
  }
  
  /**
   * Product of P(X = x|Y = x)'s induced by the given 
   * noise statistics. 
   */
  def double sumLogPrs(NoiseStatistics statistics) {
    if (!valid)
      return Double.NEGATIVE_INFINITY 
    var sum = 0.0
    sum += Math.log(p_x1_given_y0) * statistics.FN
    sum += Math.log1p(-p_x1_given_y0) * (nNegative - statistics.FN)
    sum += Math.log1p(-p_x1_given_y1) * statistics.FP
    sum += Math.log(p_x1_given_y1) * (nPositive - statistics.FP)
    return sum
  }
  
  override getCells() { binaryMatrix.cells }
  override getLoci() { binaryMatrix.loci}
  
  override logNormalization() {
    if (!valid) 
      return Double.NEGATIVE_INFINITY
    return nPositive * Math.log(1.0 - fn + fp) + nNegative * Math.log(fn + 1.0 - fp)
  }
}