package corrupt.pre

import corrupt.post.SimpleCLMatrix
import corrupt.post.CellLocusMatrix
import corrupt.GenomeMap
import corrupt.post.CLMatrixUtils
import java.io.File
import corrupt.post.ReadOnlyCLMatrix
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
import blang.inits.Implementations
import blang.inits.parsing.QualifiedName
import blang.inits.experiments.tabwriters.factories.CSV
import blang.inits.DefaultValue

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
    def  List<ReadOnlyCLMatrix> deltas(ComputeDeltas computer)
  }
  
  static val pathToSamples = "chromoplots"
  static val sampleFilePrefix = "hmms"
  static val pathToOptions = "arguments.tsv"
  static class FromPosteriorSamples implements Source {
    
    @Arg List<File> files
    
    @Arg File lociIndexFile
    
    @Arg            @DefaultValue("0.5")
    public double burnInFraction = 0.5
    
    override deltas(ComputeDeltas computer) {
      val snapshotDir = computer.results.child("snapshot")
      val indexers = loadLocusIndexers(lociIndexFile) 
      val cells = new LinkedHashSet<Cell>
      val loci = new LinkedHashSet<Locus>
      val negative = new Counter<Pair<Cell, Locus>>
      val positive = new Counter<Pair<Cell, Locus>>
      val jump = new Counter<Pair<Cell, Locus>>
      var nSamplesFull = -1
      var nSamplesPostBurnIn = -1
      for (file : files) {
        
        val arguments = CSVFile.parseTSV(new File(file, pathToOptions)).asMap()
        val cell = new Cell(arguments.get(new QualifiedName(#["cell"])).join(" "))
        cells.add(cell)
        val paths = CSV::csvFile(new File(file, pathToSamples), sampleFilePrefix)
        val states = new Counter<Locus>
        
        val sampleIndices = new LinkedHashSet<Integer>
        for (line : BriefIO::readLines(paths).indexCSV) {
          val sample = Integer.parseInt(line.get("sample"))
          sampleIndices.add(sample)
        }
        if (nSamplesFull === -1)
          nSamplesFull = sampleIndices.size
        if (nSamplesFull !== sampleIndices.size)
          throw new RuntimeException
        val cutOff = (nSamplesFull * burnInFraction) as int
        nSamplesPostBurnIn = nSamplesFull - cutOff
        val thin = Math.max(1, nSamplesFull / 5)
        
        for (line : BriefIO::readLines(paths).indexCSV) {
          val chr = Integer.parseInt(line.get("map_key_0"))
          val pos = Integer.parseInt(line.get("positions"))
          val sample = Integer.parseInt(line.get("sample"))
          val locus = indexers.get(chr).get(pos)
          loci.add(locus)
          val valueStr = line.get("value")
          val curState = if (valueStr == "NA") Double.NaN else Double.parseDouble(valueStr)
          
          if (sample % thin == 0) {
            snapshotDir.getTabularWriter("" + sample).write(
              CLMatrixUtils::CELLS -> cell,
              CLMatrixUtils::LOCI -> locus,
              CLMatrixUtils::TIP_INCL_PRS -> curState
            )
          }
          
          states.setCount(locus, curState)
          val prevLocus = indexers.get(chr).get(pos - 1)
          val prevPrevLocus = indexers.get(chr).get(pos - 2)

          if (sample > cutOff && prevLocus !== null && prevPrevLocus !== null && !Double.isNaN(curState)) {
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
        
      }
      val result = new ArrayList<ReadOnlyCLMatrix>
      for (counter : #[negative, positive, jump]) {
        val current = new SimpleCLMatrix(cells, loci)
        for (cellLocusPair : counter.keySet)
          current.set(cellLocusPair.key, cellLocusPair.value, counter.getCount(cellLocusPair) / nSamplesPostBurnIn)
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
        if (next !== null) {
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
  
  override run() {
    val matrices = source.deltas(this)
    for (var int i = 0; i < matrices.size; i++) {
      CLMatrixUtils::toCSV(matrices.get(i), results.getFileInResultFolder("matrix-" + i + ".csv.gz"))
    }
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}