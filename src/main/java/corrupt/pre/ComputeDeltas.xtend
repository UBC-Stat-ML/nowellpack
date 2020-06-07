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
import briefj.Indexer
import java.util.LinkedHashMap
import corrupt.Locus
import briefj.BriefIO
import java.util.List
import java.util.Map
import briefj.BriefMaps
import corrupt.GenomeMap.ParsedLocus
import blang.inits.experiments.Experiment
import blang.inits.Arg
import blang.inits.parsing.CSVFile
import briefj.collections.Counter
import corrupt.Cell
import java.util.LinkedHashSet
import java.util.ArrayList
import org.jgraph.graph.CellMapper
import blang.inits.Implementations
import blang.inits.parsing.QualifiedName
import blang.inits.experiments.tabwriters.factories.CSV

class ComputeDeltas extends Experiment {
  
  @Arg Source source

  def static Map<Integer, Map<Integer,Locus>> loadLocusIndexers(File indexFile) {
    val result = new LinkedHashMap<Integer, Map<Integer,Locus>>() // chr -> index -> locus
    for (line : BriefIO::readLines(indexFile).indexCSV) {
      val locus = new Locus(line.get("locus"))
      val parsed = new ParsedLocus(locus)
      val index = Integer.parseInt(line.get("indexInChr"))
      BriefMaps.getOrPutMap(result, parsed.chr).put(index, locus)
    }
    return result
  }
  
  @Implementations(FromPosteriorSamples)
  static interface Source {
    def  List<ReadOnlyCLMatrix> deltas()
  }
  
  static val pathToSamples = "results/latest/chromoplots/"
  static val sampleFilePrefix = "hmms"
  static val pathToOptions = "results/latest/arguments.tsv"
  static class FromPosteriorSamples implements Source {
    
    @Arg List<File> files
    
    @Arg File lociIndexFile
    
    override deltas() {
      val indexers = loadLocusIndexers(lociIndexFile) 
      val cells = new LinkedHashSet<Cell>
      val loci = new LinkedHashSet<Locus>
      val negative = new Counter<Pair<Cell, Locus>>
      val positive = new Counter<Pair<Cell, Locus>>
      val jump = new Counter<Pair<Cell, Locus>>
      var nSamples = -1
      for (file : files) {
        val sampleIndices = new LinkedHashSet<Integer>
        val arguments = CSVFile.parseTSV(new File(file, pathToOptions)).asMap()
        val cell = new Cell(arguments.get(new QualifiedName(#["cell"])).join(" "))
        cells.add(cell)
        val paths = CSV::csvFile(new File(file, pathToSamples), sampleFilePrefix)
        val states = new Counter<Locus>
        for (line : BriefIO::readLines(paths).indexCSV) {
          val chr = Integer.parseInt(line.get("map_key_0"))
          val pos = Integer.parseInt(line.get("positions"))
          val sample = Integer.parseInt(line.get("sample"))
          sampleIndices.add(sample)
          val locus = indexers.get(chr).get(pos)
          loci.add(locus)
          val valueStr = line.get("value")
          val curState = if (valueStr == "NA") Double.NaN else Double.parseDouble(valueStr)
          states.setCount(locus, curState)
          val prevLocus = indexers.get(chr).get(pos - 1)
          val prevPrevLocus = indexers.get(chr).get(pos - 2)
          if (prevLocus !== null && prevPrevLocus !== null && !Double.isNaN(curState)) {
            val s0 = states.getCount(prevPrevLocus)
            val s1 = states.getCount(prevLocus)
            val s2 = curState
            val delta = 
              if (!Double.isNaN(s1)) {
                // normal jump
                s2 - s1
              } else { // => s0 NA s1
                // check
                if (Double.isNaN(s0)) throw new RuntimeException
                s2 - s0
              }
            if (delta < 0)
                negative.incrementCount(cell -> locus, 1.0)
            if (delta > 0)
              positive.incrementCount(cell -> locus, 1.0)
            if (delta == 0 && Double.isNaN(s1))
              jump.incrementCount(cell -> locus, 1.0)
          }
        }
        if (nSamples === -1)
          nSamples = sampleIndices.size
        if (nSamples !== sampleIndices.size)
          throw new RuntimeException
      }
      val result = new ArrayList<ReadOnlyCLMatrix>
      for (counter : #[negative, positive, jump]) {
        val current = new SimpleCLMatrix(cells, loci)
        for (cellLocusPair : counter.keySet)
          current.set(cellLocusPair.key, cellLocusPair.value, counter.getCount(cellLocusPair) / nSamples)
        result.add(ReadOnlyCLMatrix::readOnly(current))
      }
      return result
    }
    
  }
  
  /**
   * For a single matrix, e.g. output of HMMCopy
   */
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
  
//  def static void main(String [] args) {
//    val cns = CLMatrixUtils::fromCSV(new File(
//      //"/Users/bouchard/experiments/sitka-nextflow/data/raw/535/cn.csv.gz"))
//      "/Users/bouchard/experiments/corrupt-nextflow/analysis-chromobreak/data/SA535_uncor_state.csv"))
//    
//    val deltas = delta(cns)
//    val neg = ReadOnlyCLMatrix::readOnly(deltas.key)
//    val pos = ReadOnlyCLMatrix::readOnly(deltas.value)
//    val phylo = PerfectPhylo::parseNewick(new File("/Users/bouchard/experiments/corrupt-nextflow/results/all/2020-06-04-08-23-49-srdge3qP.exec/consensus.newick"))
////    val viz = new PerfectPhyloViz(phylo, #[neg, pos], Viz::fixHeight(300))
////    viz.output(new File("temp2.pdf"))
//    val dir = new File("temp3")
//    PerfectPhyloViz::visualizePerChromosome(dir, phylo, #[neg, pos], Viz::fixHeight(300))
//  }
  
  override run() {
    val matrices = source.deltas
    for (var int i = 0; i < matrices.size; i++) {
      CLMatrixUtils::toCSV(matrices.get(i), results.getFileInResultFolder("matrix-" + i + ".csv.gz"))
    }
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}