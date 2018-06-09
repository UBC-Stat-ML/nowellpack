package corrupt.post

import com.rits.cloning.Immutable
import org.eclipse.xtend.lib.annotations.Data
import corrupt.Cell
import corrupt.Locus
import org.eclipse.xtend.lib.annotations.Delegate
import blang.inits.DesignatedConstructor
import blang.types.Plate
import blang.types.Plated
import blang.inits.ConstructorArg

@Immutable
@Data class ReadOnlyCLMatrix implements CellLocusMatrix {
  @Delegate
  val CellLocusMatrix enclosed
  
  @DesignatedConstructor
  public static def ReadOnlyCLMatrix create(
      @ConstructorArg("cells") Plate<Cell> cells, 
      @ConstructorArg("loci") Plate<Locus> loci, 
      @ConstructorArg("tipInclusionProbabilities") Plated<Double> inclusionPrs
  ) { 
    println("Loading tip inclusion probabilities")
    val enclosed = new SimpleCLMatrix(cells.indices.map[key].toSet, loci.indices.map[key].toSet)
    for (cell : cells.indices) 
      for (locus : loci.indices) 
        enclosed.setTip(cell.key, locus.key, inclusionPrs.get(cell, locus))
    return new ReadOnlyCLMatrix(enclosed)   
  } 
}