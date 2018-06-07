package corrupt.post

class CLMatrixUtils {
  def static double distance(CellLocusMatrix cl1, CellLocusMatrix cl2) {
    checkCompatible(cl1, cl2)
    var sum = 0.0
    for (cell : cl1.cells)
      for (locus : cl1.loci)
        sum += Math.abs(cl1.getTipAsDouble(cell, locus) - cl2.getTipAsDouble(cell, locus))
    return sum / (cl1.cells.size * cl2.cells.size)
  }
  def static void checkCompatible(CellLocusMatrix cl1, CellLocusMatrix cl2) {
    if (cl1.cells != cl2.cells || cl1.loci != cl2.loci)
      throw new RuntimeException
  }
}