package corrupt.pre

import corrupt.post.SimpleCLMatrix
import corrupt.post.CellLocusMatrix
import corrupt.GenomeMap
import corrupt.post.CLMatrixUtils
import java.io.File
import corrupt.PerfectPhylo
import corrupt.viz.PerfectPhyloViz
import viz.core.Viz
import corrupt.post.ReadOnlyCLMatrix

class ComputeDeltas {
  
  def static Pair<SimpleCLMatrix,SimpleCLMatrix> delta(CellLocusMatrix matrix) {
    val negative = new SimpleCLMatrix(matrix.cells, matrix.loci)
    val positive = new SimpleCLMatrix(matrix.cells, matrix.loci)
    val map = new GenomeMap(matrix.loci)
    for (chr : map.orderedChromosomes) {
      val orderedLoci = map.orderedLoci(chr)
      for (var int i = 0; i < orderedLoci.size - 1; i++) {
        val locus = orderedLoci.get(i)
        val next = orderedLoci.get(i + 1)
        if (next !== null && GenomeMap.isAdjacent(locus, next)) {
          for (cell : matrix.cells) {
            val delta = matrix.get(cell, next) - matrix.get(cell, locus)
            if (delta < 0)
              negative.set(cell, locus, -delta)
            if (delta > 0)
              positive.set(cell, locus, delta)
          }
        }
      }
    }
    return negative -> positive
  }
  
  def static void main(String [] args) {
    val cns = CLMatrixUtils::fromCSV(new File(
      //"/Users/bouchard/experiments/sitka-nextflow/data/raw/535/cn.csv.gz"))
      "/Users/bouchard/experiments/corrupt-nextflow/analysis-chromobreak/data/SA535_uncor_state.csv"))
    
    val deltas = delta(cns)
    val neg = ReadOnlyCLMatrix::readOnly(deltas.key)
    val pos = ReadOnlyCLMatrix::readOnly(deltas.value)
    val phylo = PerfectPhylo::parseNewick(new File("/Users/bouchard/experiments/corrupt-nextflow/results/all/2020-06-04-08-23-49-srdge3qP.exec/consensus.newick"))
//    val viz = new PerfectPhyloViz(phylo, #[neg, pos], Viz::fixHeight(300))
//    viz.output(new File("temp2.pdf"))
    val dir = new File("temp3")
    PerfectPhyloViz::visualizePerChromosome(dir, phylo, #[neg, pos], Viz::fixHeight(300))
  }
}