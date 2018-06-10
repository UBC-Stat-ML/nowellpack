package corrupt.post

import blang.inits.DesignatedConstructor
import blang.inits.Input
import com.rits.cloning.Immutable
import java.io.File
import org.eclipse.xtend.lib.annotations.Data
import org.eclipse.xtend.lib.annotations.Delegate

@Immutable
@Data class ReadOnlyCLMatrix implements CellLocusMatrix {
  @Delegate
  val CellLocusMatrix enclosed
  
  @DesignatedConstructor
  public static def ReadOnlyCLMatrix create(
      @Input() String path
  ) { 
    println("Loading tip inclusion probabilities")
    return new ReadOnlyCLMatrix(CLMatrixUtils::fromCSV(new File(path)))    
  } 
}