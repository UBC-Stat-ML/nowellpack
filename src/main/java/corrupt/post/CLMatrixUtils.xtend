package corrupt.post

import bayonet.distributions.Random
import corrupt.PerfectPhylo
import blang.distributions.Normal

import static blang.types.StaticUtils.*
import bayonet.distributions.Multinomial
import briefj.BriefIO
import java.io.File
import java.util.ArrayList
import corrupt.Locus
import corrupt.Cell
import blang.inits.providers.CoreProviders

class CLMatrixUtils {
  
  public static val String CELLS = "cells"
  public static val String LOCI = "loci"
  public static val String TIP_INCL_PRS = "tipInclusionProbabilities"
  
  static def void toCSV(CellLocusMatrix matrix, File f) {
    val out = BriefIO.output(f)
    out.println('''«CELLS»,«LOCI»,«TIP_INCL_PRS»''')
    for (cell : matrix.cells)
      for (locus : matrix.loci)
        out.println('''«cell»,«locus»,«matrix.getTipAsDouble(cell, locus)»''')
    out.close
  }
  
  static def SimpleCLMatrix fromCSV(File f) {
    val cells = new ArrayList<Cell>
    val loci = new ArrayList<Locus>
    val values = new ArrayList<Double>
    for (line : BriefIO::readLines(f).indexCSV) {
      if (!line.containsKey(CELLS) || !line.containsKey(LOCI) || !line.containsKey(TIP_INCL_PRS))
        throw new RuntimeException
      cells.add(new Cell(line.get(CELLS)))
      loci.add(new Locus(line.get(LOCI)))
      values.add(CoreProviders::parse_double(line.get(TIP_INCL_PRS)))
    }
    val result = new SimpleCLMatrix(cells, loci)
    for (var int i = 0; i < cells.size; i++) 
      result.setTip(cells.get(i), loci.get(i), values.get(i))
    return result
  }
  
  static def ReadOnlyCLMatrix syntheticInclusionPrs(Random rand, PerfectPhylo phylo, double stdDev) {
    val syntheticInclusionPrs = new SimpleCLMatrix(phylo.cells, phylo.loci)
    val inclNormal = Normal::distribution(fixedReal(0.0), fixedReal(stdDev * stdDev))
    val exclNormal = Normal::distribution(fixedReal(1.0), fixedReal(stdDev * stdDev))
    for (locus : phylo.loci) {
      val tips = phylo.getTips(locus)
      for (cell : phylo.cells) {
        val dist = if (tips.get(cell)) inclNormal else exclNormal
        val observation = dist.sample(rand)
        val prs = newDoubleArrayOfSize(2)
        prs.set(0, exclNormal.logDensity(observation))
        prs.set(1, inclNormal.logDensity(observation))
        Multinomial::expNormalize(prs)
        syntheticInclusionPrs.setTip(cell, locus, prs.get(1)) 
      }
    }
    return new ReadOnlyCLMatrix(syntheticInclusionPrs)
  }
  
  static def void checkCompatible(CellLocusMatrix cl1, CellLocusMatrix cl2) {
    if (cl1.cells != cl2.cells || cl1.loci != cl2.loci)
      throw new RuntimeException
  }
}