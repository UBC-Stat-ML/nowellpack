package chromobreak

import java.io.File
import blang.inits.experiments.Experiment
import blang.inits.Arg
import blang.inits.DefaultValue
import corrupt.post.CLMatrixUtils
import corrupt.GenomeMap

/**
 * This will:
 * - create tidy formats
 * - re-index loci (and cell if requested)
 * - create index files
 * - can split into cell specific directories too
 */
class TidifyCounts extends Experiment {
  
  @Arg
  public File countFile
  
  @Arg @DefaultValue("INF")
  public int nCells = Integer.MAX_VALUE
  
  @Arg @DefaultValue("INF")
  public int nChromosomes = Integer.MAX_VALUE
  
  @Arg(description = "Use integer indices for cells") @DefaultValue("false")
  public boolean useInteger = false // set to true for legacy
  
  @Arg                @DefaultValue("false")
  public boolean cellSpecificFiles = false
  
  public static String DATA_OUTPUT = "tidy"
  
  override run() {
    val parsed = CLMatrixUtils::fromCSV(countFile)
    val map = new GenomeMap(parsed.loci)
    
    if (!map.lociAdjacent)
      throw new RuntimeException("We assume all loci in a chromosome are adjacent after proper sorting. Are you running on a subset of the loci?")
    
    val cellList = parsed.cells.toList    
    for (var int cellIndex = 0; cellIndex < cellList.size; cellIndex++) 
      if (cellIndex <= nCells) {
        val cell = cellList.get(cellIndex)
        for (chr : map.orderedChromosomes) 
          if (chr < nChromosomes) {
            val lociList = map.orderedLoci(chr)
            for (var int posIndex = 0; posIndex < lociList.size; posIndex++)  {
              val locus = lociList.get(posIndex)
              // re-index by position
              val Object cellCode = if (useInteger) cellIndex else cell
              val writer = 
                if (cellSpecificFiles) {
                  results.child(DATA_OUTPUT).child(cellCode.toString).getTabularWriter("data")
                } else 
                  results.getTabularWriter(DATA_OUTPUT).child("cell", cellCode)
              
              writer.write(
                "chromosomes" -> chr,
                "positions" -> posIndex,
                "value" -> parsed.get(cell, locus)
              )
              
              // write the index for one cell only
              if (cellIndex === 0) {
                results.getTabularWriter("lociIndex").write(
                  "locus" -> locus,
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