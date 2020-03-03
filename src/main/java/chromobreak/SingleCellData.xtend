package chromobreak

import blang.types.Plate
import blang.io.GlobalDataSource
import blang.inits.Arg
import blang.core.IntVar
import blang.types.Plated
import blang.core.RealVar
import blang.runtime.internals.objectgraph.SkipDependency
import blang.inits.DesignatedConstructor

class SingleCellData {
  @Arg public GlobalDataSource source
  
  @Arg public Plate<String> chromosomes
  
  @SkipDependency(isMutable = false)
  @Arg public Plate<Integer> positions
  
  @Arg public Plated<RealVar> gcContents
  @Arg public Plated<IntVar> readCounts
  
  @DesignatedConstructor
  new () {
    ChromoPostProcessor.data = this
  }
}