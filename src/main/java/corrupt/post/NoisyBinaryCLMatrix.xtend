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
  
  @DesignatedConstructor
  new(
    @ConstructorArg("binaryMatrix") ReadOnlyCLMatrix binaryMatrix, 
    @ConstructorArg("fpRate") double fpRate, 
    @ConstructorArg("fnRate") double fnRate
  ) {
    this.binaryMatrix = binaryMatrix
    this.fpRate = new RealConstant(fpRate)
    this.fnRate = new RealConstant(fnRate)
  }
  
  new(
    ReadOnlyCLMatrix binaryMatrix, 
    RealVar fpRate, 
    RealVar fnRate
  ) {
    this.binaryMatrix = binaryMatrix
    this.fpRate = fpRate
    this.fnRate = fnRate 
  }
  
  override getTipAsDouble(Cell cell, Locus locus) { 
    val boolValue = 
      switch(binaryMatrix.getTipAsDouble(cell, locus)) {
        case 0.0 : false
        case 1.0 : true
        default  : throw new RuntimeException
      }
    val fp = fpRate.doubleValue
    val fn = fnRate.doubleValue
    if (fp < 0.0 || fp > 1.0 || fn < 0.0 || fn > 1.0)
      throw new RuntimeException
    if (boolValue) return (1.0 - fn) / (1.0 - fn + fp)
    else           return fn / (fn + 1.0 - fp)
  }
  
  override getCells() { binaryMatrix.cells }
  override getLoci() { binaryMatrix.loci}
}