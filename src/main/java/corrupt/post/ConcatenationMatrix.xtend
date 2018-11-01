package corrupt.post

import corrupt.Cell
import corrupt.Locus
import java.util.List
import org.eclipse.xtend.lib.annotations.Data
import java.util.Map
import java.util.HashMap
import java.util.LinkedHashSet

@Data class ConcatenationMatrix implements CellLocusMatrix {
  
  val List<CellLocusMatrix> matrices
  val Map<Locus,Integer> index = new HashMap
  val union = new LinkedHashSet<Locus>
  
  new (List<CellLocusMatrix> matrices) {
    if (matrices.empty)
      throw new RuntimeException
    this.matrices = matrices
    // check cell sets match
    for (var int i = 1; i < matrices.size; i++)
      if (matrices.get(0).cells != matrices.get(1).cells)
        throw new RuntimeException("The set of cells need to match when performing concatenation analyses.")
    // check loci disjoint
    
    var sum = 0
    for (matrix : matrices) {
      union.addAll(matrix.loci)
      sum += matrix.loci.size
    }
    if (sum != union.size)
      throw new RuntimeException("Loci in any two concatenated matrices should bear distinct identifiers.")
    // create index
    for (var int i = 0; i < matrices.size; i++)
      for (locus : matrices.get(i).loci)
        index.put(locus, i)
  }
  
  override get(Cell cell, Locus locus) {
    matrices.get(index.get(locus)).get(cell, locus)
  }
  
  override getCells() {
    matrices.get(0).cells
  }
  
  override getLoci() {
    union
  }
  
}