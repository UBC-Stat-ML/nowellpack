package corrupt.post

import corrupt.Locus
import corrupt.Cell
import java.util.Map

class NoiseStatistics {
    public var nOnes  = 0.0  // # latent indicators with value one
    public var nZeros = 0.0  // ... with value zero
    public var FP = 0.0
    public var FN = 0.0
    def fpRate() { FP / nZeros }
    def fnRate() { FN / nOnes }
    def void add(Locus locus, Map<Cell, Boolean> tips, CellLocusMatrix binaryMatrix)      { add(locus, tips, binaryMatrix, 1)}
    def void subtract(Locus locus, Map<Cell, Boolean> tips, CellLocusMatrix binaryMatrix) { add(locus, tips, binaryMatrix, -1)}
    def void add(Locus locus, Map<Cell, Boolean> tips, CellLocusMatrix binaryMatrix, int increment) {
      for (entry : tips.entrySet) {
        val truth = entry.value
        val noisy = 
          switch (binaryMatrix.get(entry.key, locus)) {
            case 1.0 : true
            case 0.0 : false
            default  : throw new RuntimeException
          }
        if (truth && !noisy) FN += increment
        if (!truth && noisy) FP += increment
        if (truth) nOnes += increment
        else       nZeros += increment
      }
    }
  }