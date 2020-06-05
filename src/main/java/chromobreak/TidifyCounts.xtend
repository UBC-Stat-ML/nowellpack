package chromobreak

import java.io.File
import blang.inits.experiments.Experiment
import blang.inits.Arg
import blang.inits.DefaultValue
import corrupt.post.CLMatrixUtils
import corrupt.GenomeMap
import corrupt.GenomeMap.ParsedLocus

/**
 * This will also
 * - re-index loci (and cell if requested)
 */
class TidifyCounts extends Experiment {
  
  @Arg
  public File countFile
  
  @Arg @DefaultValue("INF")
  public int nCells = Integer.MAX_VALUE
  
  @Arg(description = "Use integer indices for cells") @DefaultValue("false")
  public boolean useInteger = false // set to true for legacy
  
  override run() {
    val parsed = CLMatrixUtils::fromCSV(countFile)
    val map = new GenomeMap(parsed.loci)
    
    val cellList = parsed.cells.toList    
    for (var int cellIndex = 0; cellIndex < cellList.size; cellIndex++) 
      if (cellIndex < nCells) {
        val cell = cellList.get(cellIndex)
        for (chr : map.orderedChromosomes) {
          val lociList = map.orderedLoci(chr)
          for (var int posIndex = 0; posIndex < lociList.size; posIndex++)  {
            val locus = lociList.get(posIndex)
            // re-index by position
            results.getTabularWriter("tidy").write(
              "cell" -> if (useInteger) cellIndex else cell,
              "chromosomes" -> chr,
              "positions" -> posIndex,
              "value" -> parsed.get(cell, locus)
            )
            
            // write the index for one cell only
            if (cellIndex === 0) {
              val parsedLocus = new ParsedLocus(locus)
              results.getTabularWriter("lociIndex").write(
                "chromosomes" -> chr,
                "start" -> parsedLocus.leftOneIndexedIncl,
                "end" -> parsedLocus.rightOneIndexedIncl,
                "indexInChr" -> posIndex
              )
            }
            
            if (useInteger) {
              results.getTabularWriter("cellsIndex").write(
                "cell" -> cell,
                "integerId" -> cellIndex
              )
            }
          }
        }
    }
  }
  
  def static void main(String [] args){
    Experiment.startAutoExit(args)
  }

}