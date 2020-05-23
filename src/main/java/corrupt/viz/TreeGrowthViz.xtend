package corrupt.viz

import blang.inits.experiments.Experiment
import briefj.BriefIO
import corrupt.PerfectPhylo
import corrupt.post.ReadOnlyCLMatrix
import java.util.Optional
import corrupt.post.CellLocusMatrix
import corrupt.post.SimpleCLMatrix
import blang.inits.Arg
import java.io.File
import java.util.List
import viz.core.PublicSize
import java.util.ArrayList
import blang.inits.DefaultValue

class TreeGrowthViz extends Experiment { 
  @Arg(description = "CSV of phylogenies") File phylogenies 
  @Arg List<ReadOnlyCLMatrix> matrices 
  @Arg PublicSize size
  @Arg Optional<PerfectPhylo> refPhylo = Optional.empty
  @Arg Optional<List<Integer>> codes = Optional.empty 
  
  @Arg      @DefaultValue("value")
  String phylogeniesColName = "value"
  
  @Arg        @DefaultValue("iteration")
  String iterationColName = "iteration"
  
  override run() {
    var List<ReadOnlyCLMatrix> prevMatrices = null
    for (line : BriefIO::readLines(phylogenies).indexCSV) {
      val tree = new PerfectPhylo(line.get(phylogeniesColName))
      val List<ReadOnlyCLMatrix> restricedMatrices = new ArrayList
      for (m : matrices)
        restricedMatrices.add(ReadOnlyCLMatrix::readOnly(eraseUnassigned(m, tree)))
      new PerfectPhyloViz(
        tree, restricedMatrices, size, refPhylo, codes
      ).output(results.getFileInResultFolder(line.get(iterationColName) + ".pdf"))
      
      if (prevMatrices !== null) {
        val List<ReadOnlyCLMatrix> deltas = new ArrayList
        for (var int i = 0; i < prevMatrices.size; i++) 
          deltas.add(ReadOnlyCLMatrix::readOnly(delta(restricedMatrices.get(i), prevMatrices.get(i))))
        new PerfectPhyloViz(
          tree, deltas, size, refPhylo, codes
        ).output(results.getFileInResultFolder("delta_" + line.get(iterationColName) + ".pdf"))
      }
      
      prevMatrices = restricedMatrices
    }
  }
  
  def static CellLocusMatrix eraseUnassigned(CellLocusMatrix m, PerfectPhylo phylo) {
    val result = new SimpleCLMatrix(m.cells, m.loci)
    for (locus : m.loci) {
      val current = phylo.getTips(locus)
      var found = false
      for (value : current.values)
        if (value)
          found = true
      if (found) 
        for (cell : m.cells)
          result.set(cell, locus, m.get(cell, locus))
    }
    return result
  }
  
  def static CellLocusMatrix delta(CellLocusMatrix current, CellLocusMatrix prev) {
    val result = new SimpleCLMatrix(current.cells, current.loci) 
    for (locus : current.loci)
      for (cell : current.cells)
        result.set(cell, locus, Math.abs(current.get(cell, locus) - prev.get(cell, locus)))
    return result
  }
  
  static def void main(String [] args) { 
    Experiment.startAutoExit(args)
  }
}