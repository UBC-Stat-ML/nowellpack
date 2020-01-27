package corrupt

import briefj.BriefIO
import java.io.File

class Old {
  def static void main(String [] args){
    val nCells = 500
    val nChr = 20
    val outPath = "/Users/bouchard/experiments/corrupt-nextflow/analysis-raw-data/temp.csv"
    val out = BriefIO.output(new File(outPath))
    out.println("cell,chr,pos,norm_reads")
    for (line : BriefIO.readLines("/Users/bouchard/experiments/corrupt-nextflow/analysis-raw-data/SA535_cn_reads.csv").splitCSV.skip(1)) {
      val chr = GenomeMap::chromosomeIndex(line.get(0))
      if (chr <= nChr)
        for (c : 0 ..< nCells)
          out.println("cell_" + c + "," + chr + "," + line.get(1) + "," + line.get(3+c))
    }
    out.close
  }

}