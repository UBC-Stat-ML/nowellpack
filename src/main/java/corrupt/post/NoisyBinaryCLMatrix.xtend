package corrupt.post

import org.eclipse.xtend.lib.annotations.Data
import corrupt.Cell
import corrupt.Locus
import blang.core.RealVar
import blang.inits.ConstructorArg
import blang.inits.DesignatedConstructor
import xlinear.Matrix
import java.util.List
import java.util.Map
import static xlinear.MatrixOperations.*
import static extension xlinear.MatrixExtensions.*
import java.util.HashMap
import blang.core.RealConstant
import blang.core.WritableRealVar
import blang.runtime.internals.objectgraph.SkipDependency
import xlinear.internals.CommonsDenseMatrix

@Data class NoisyBinaryCLMatrix implements CellLocusMatrix {
  val BinaryCLMatrix binaryMatrix 
  val Matrix fpRates
  val Matrix fnRates
  @SkipDependency(isMutable = false) val Matrix nPositives
  @SkipDependency(isMutable = false) val Matrix nNegatives
  val boolean global
  @SkipDependency(isMutable = false) val Map<Locus, Integer> parametersMap
  
  def int parameter(Locus locus) {
    if (global)
      return 0
    else
      return parametersMap.get(locus)
  }
  
  @DesignatedConstructor
  new(
    @ConstructorArg("binaryMatrix") BinaryCLMatrix binaryMatrix, 
    @ConstructorArg("fpRate") double fpRate, 
    @ConstructorArg("fnRate") double fnRate
  ) {
    this(binaryMatrix, new RealConstant(fpRate), new RealConstant(fnRate)) 
  }
  
  new(
    BinaryCLMatrix binaryMatrix, 
    RealVar fpRate, 
    RealVar fnRate
  ) {
    this(binaryMatrix, viewScalarAsMatrix(fpRate), viewScalarAsMatrix(fnRate))
  }

  new(
    BinaryCLMatrix binaryMatrix, 
    Matrix fpRates, 
    Matrix fnRates
  ) {
    if (fpRates.nEntries !== fnRates.nEntries || !fpRates.vector || !fnRates.vector)
      throw new RuntimeException
    this.global = fpRates.nEntries === 1
    this.binaryMatrix = binaryMatrix
    this.fpRates = fpRates
    this.fnRates = fnRates
    this.parametersMap = new HashMap 
    this.nPositives = dense(nParams)
    this.nNegatives = dense(nParams)
    var int index = 0
    for (locus : binaryMatrix.loci) {
      if (!global)
        parametersMap.put(locus, index)
      for (cell : binaryMatrix.cells)
        if (binaryValue(cell, locus))
          nPositives.increment(index, 1)
        else
          nNegatives.increment(index, 1)
      if (!global)
        index++
    }
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
    if (!valid) 
      return Double::NEGATIVE_INFINITY 
    val boolValue = binaryValue(cell, locus)
    val paramIndex = parameter(locus)
    if (boolValue) return p_x1_given_y1(paramIndex)
    else           return p_x1_given_y0(paramIndex)
  }
  
  def boolean valid() {
    for (p : 0 ..< nParams)   
      if (fp(p) <= 0.0 || fp(p) >= 1.0 || fn(p) <= 0.0 || fn(p) >= 1.0)
        return false
    return true
  }
  
  def double fp(int paramIndex) { fpRates.get(paramIndex) }
  def double fn(int paramIndex) { fnRates.get(paramIndex) }
  
  def double nNegative(int paramIndex) { nNegatives.get(paramIndex) }
  def double nPositive(int paramIndex) { nPositives.get(paramIndex) }
  
  /**
   * P(X = 1 | Y = 1)
   */
  def double p_x1_given_y1(int paramIndex) {
    return (1.0 - fn(paramIndex)) / (1.0 - fn(paramIndex) + fp(paramIndex))
  }
  
  /**
   * P(X = 1 | Y = 0)
   */
  def double p_x1_given_y0(int paramIndex) {
    return fn(paramIndex) / (fn(paramIndex) + 1.0 - fp(paramIndex))
  }
  
  /**
   * Product of P(X = x|Y = x)'s induced by the given 
   * noise statistics. 
   */
  def double sumLogPrs(List<NoiseStatistics> statistics) {
    if (!valid)
      return Double.NEGATIVE_INFINITY 
    if (statistics.size != nParams)
      throw new RuntimeException
    var sum = 0.0
    for (p : 0 ..< nParams) {
      sum += Math.log(p_x1_given_y0(p)) * statistics.get(p).FN
      sum += Math.log1p(-p_x1_given_y0(p)) * (nNegative(p) - statistics.get(p).FN)
      sum += Math.log1p(-p_x1_given_y1(p)) * statistics.get(p).FP
      sum += Math.log(p_x1_given_y1(p)) * (nPositive(p) - statistics.get(p).FP)
    }
    return sum
  }
  
  override getCells() { binaryMatrix.cells }
  override getLoci() { binaryMatrix.loci}
  
  def int nParams() { if (global) 1 else loci.size }
  
  override logNormalization() {
    if (!valid) 
      return Double.NEGATIVE_INFINITY
    var sum = 0.0
    for (p : 0 ..< nParams)
      sum += nPositive(p) * Math.log(1.0 - fn(p) + fp(p)) + nNegative(p) * Math.log(fn(p) + 1.0 - fp(p))
    return sum
  }
    
  private def static Matrix viewScalarAsMatrix(RealVar rv) {
    val result = new CommonsDenseMatrix(null) {
      override nRows() { 1 }
      override nCols() { 1 }
      override get(int row, int col) { rv.doubleValue } 
      override set(int row, int col, double v) { (rv as WritableRealVar).set(v) }
    }
    if (rv instanceof WritableRealVar)
      return result
    else
      return result.readOnlyView
  }
}