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
import xlinear.Matrix
import blang.distributions.Generators
import java.util.HashSet
import corrupt.distances.LocusEdgeStat
import java.util.Set
import java.util.LinkedHashSet
import java.util.List
import corrupt.GenomeMap

class CLMatrixUtils {
  
  public static val String CELLS = "cells"
  public static val String LOCI = "loci"
  public static val String TIP_INCL_PRS = "tipInclusionProbabilities"
  
  static def void toCSV(CellLocusMatrix matrix, File f) {
    toCSV(matrix, f, null)
  }
  static def void toCSV(CellLocusMatrix matrix, File f, CellLocusMatrix reference) {
    val out = BriefIO.writer(f)
    out.append('''«CELLS»,«LOCI»,«TIP_INCL_PRS»«IF reference !== null»,reference«ENDIF»''' + "\n")
    for (cell : matrix.cells)
      for (locus : matrix.loci)
        out.append('''«cell»,«locus»,«matrix.get(cell, locus)»«IF reference !== null»,«reference.get(cell, locus)»«ENDIF»''' + "\n")
    out.close
  }
  
  static def int locusEdgeDistance(PerfectPhylo phylo1, PerfectPhylo phylo2) {
    val s1 = LocusEdgeStat::locusEdges(phylo1.collapsedTree)
    val s2 = LocusEdgeStat::locusEdges(phylo2.collapsedTree)
    val union = new HashSet => [
      addAll(s1)
      addAll(s2)
    ]
    val intersection = new HashSet(s1)
    intersection.retainAll(s2)
    union.removeAll(intersection)
    return union.size
  }
  
  static def double distance(PerfectPhylo phylo1, PerfectPhylo phylo2) {
    return distance(fromPhylo(phylo1), fromPhylo(phylo2))
  }
  
  static def double distance(SimpleCLMatrix mtx1, SimpleCLMatrix mtx2) {
    checkCompatible(mtx1, mtx2)
    return distance(mtx1.matrix, mtx2.matrix)
  }
  
  static def double distance(Matrix mtx1, Matrix mtx2) {
    return delta(mtx1, mtx2) / mtx1.nEntries
  }
  
  static def double delta(Matrix mtx1, Matrix mtx2) {
    val diff = mtx1 - mtx2
    return diff.nonZeroEntries().map[Math.abs(it)].sum()
  }
  
  static def SimpleCLMatrix fromCSV(File f) {
    if (isLegacyFormat(headers(f)))
      return fromLegacyCSV(f)
    else 
      return fromTidyCSV(f)
  }
  
  static private def List<String> headers(File f) {
    BriefIO.readLines(f).splitCSV.first.get.toList
  }
  
  /**
   * Assumes 
   * chr,start,end,[name of first cell], ...
   */
  static def SimpleCLMatrix fromLegacyCSV(File f) {
    val headers = headers(f)
    if (!isLegacyFormat(headers)) throw new RuntimeException
    val cells = headers.subList(3, headers.size).map[new Cell(it.trim)].toList
    if (cells.size !== cells.toSet.size) 
      throw new RuntimeException("Duplicate cell: " + cells)
    val loci = new LinkedHashSet<Locus>
    for (line : BriefIO.readLines(f).splitCSV.skip(1)) {
      val locus = parseLocus(line)
      if (loci.contains(locus)) throw new RuntimeException("Duplicate locus: " + locus)
      loci.add(locus)
    }
    val result = new SimpleCLMatrix(cells, loci)
    for (line : BriefIO.readLines(f).splitCSV.skip(1)) {
      val locus = parseLocus(line)
      for (var int i = 0; i < cells.size; i++) {
        val cell = cells.get(i)
        val value = Double::parseDouble(line.get(i+3))
        result.set(cell, locus, value)
      }
    }
    return result
  }
  
  def static Locus parseLocus(List<String> line) {
    val chr = line.get(0)
    val left = Integer::parseInt(line.get(1))
    val right = Integer::parseInt(line.get(2))
    return GenomeMap::locus(chr, left, right)
  }
  
  static private def boolean isLegacyFormat(List<String> headers) {
    if (headers.size < 3) return false
    return headers.subList(0, 3).equals(#["chr","start","end"])
  }
  
  /**
   * Assumes tidy format where one 
   */
  static def SimpleCLMatrix fromTidyCSV(File f) {
    val fields = BriefIO.readLines(f).indexCSV.first.get
    val cellsField = pickExactlyOne(fields.keySet, CELLS, "cell")
    val lociField = pickExactlyOne(fields.keySet, LOCI, "locus")
    val tipIncl = pickExactlyOne(fields.keySet, TIP_INCL_PRS, "value")
    return fromTidyCSV(f, cellsField, lociField, tipIncl)
  }
  
  static private def pickExactlyOne(Set<String> fields, String candidate1, String candidate2) {
    val c1 = fields.contains(candidate1)
    val c2 = fields.contains(candidate2)
    if (c1 == c2) throw new RuntimeException("Should contain fields matching exactly " + candidate1 + " or " + candidate2)
    if (c1) candidate1 else candidate2
  }
  
  static def SimpleCLMatrix fromTidyCSV(File f, String CELLS, String LOCI, String TIP_INCL_PRS) {
    val cells = new ArrayList<Cell>
    val loci = new ArrayList<Locus>
    val values = new ArrayList<Double>
    for (line : BriefIO::readLines(f).indexCSV) {
      if (!line.containsKey(CELLS) || !line.containsKey(LOCI) || !line.containsKey(TIP_INCL_PRS))
        throw new RuntimeException("The following column headers should be present in " + f.name + ": " + 
          CELLS + ", " + LOCI + ", " + TIP_INCL_PRS + "\nFound: " + line.keySet.join(", ")
        )
      cells.add(new Cell(line.get(CELLS)))
      loci.add(new Locus(line.get(LOCI)))
      values.add(CoreProviders::parse_double(line.get(TIP_INCL_PRS)))
    }
    val result = new SimpleCLMatrix(cells, loci)
    for (var int i = 0; i < cells.size; i++) { 
      val cell = cells.get(i)
      val locus = loci.get(i)
      if (result.get(cell, locus) !== 0.0)
        throw new RuntimeException("Duplicate entry for cell=" + cell + ",locus=" + locus + " in " + f.absolutePath)
      result.set(cell, locus, values.get(i))
    }
    return result
  }
  
  static def SimpleCLMatrix fromPhylo(PerfectPhylo phylo) {
    val result = new SimpleCLMatrix(phylo.cells, phylo.loci)
    result += phylo 
    return result
  }
  
  static def ReadOnlyCLMatrix syntheticInclusionPrs(Random rand, PerfectPhylo phylo, double stdDev) {
    syntheticInclusionPrs(rand, phylo, stdDev, phylo.loci)
  }
  static def ReadOnlyCLMatrix syntheticInclusionPrs(Random rand, PerfectPhylo phylo, double stdDev, Set<Locus> loci) {
    val syntheticInclusionPrs = new SimpleCLMatrix(phylo.cells, loci)
    val inclNormal = Normal::distribution(fixedReal(0.0), fixedReal(stdDev * stdDev))
    val exclNormal = Normal::distribution(fixedReal(1.0), fixedReal(stdDev * stdDev))
    for (locus : loci) {
      val tips = phylo.getTips(locus)
      for (cell : phylo.cells) {
        val dist = if (tips.get(cell)) inclNormal else exclNormal
        val observation = dist.sample(rand)
        val prs = newDoubleArrayOfSize(2)
        prs.set(0, exclNormal.logDensity(observation))
        prs.set(1, inclNormal.logDensity(observation))
        Multinomial::expNormalize(prs)
        syntheticInclusionPrs.set(cell, locus, prs.get(1)) 
      }
    }
    return ReadOnlyCLMatrix.readOnly(syntheticInclusionPrs)
  }
  
  static def BinaryCLMatrix syntheticPerturbedBinaryMatrix(Random rand, PerfectPhylo phylo, double fpRate, double fnRate) {
    syntheticPerturbedBinaryMatrix(rand, phylo, fpRate, fnRate, phylo.loci)
  }
  static def BinaryCLMatrix syntheticPerturbedBinaryMatrix(Random rand, PerfectPhylo phylo, double fpRate, double fnRate, Set<Locus> loci) {
    val syntheticInclusionPrs = new SimpleCLMatrix(phylo.cells, loci)
    for (locus : loci) {
      val tips = phylo.getTips(locus)
      for (cell : phylo.cells) {
        var indic = tips.get(cell)
        if (indic) {
          if (Generators::bernoulli(rand, fnRate))
            indic = false
        } else {
          if (Generators::bernoulli(rand, fpRate))
            indic = true
        }
        syntheticInclusionPrs.set(cell, locus, if (indic) 1.0 else 0.0) 
      }
    }
    return BinaryCLMatrix::binary(syntheticInclusionPrs)
  }
  
  static def void checkCompatible(SimpleCLMatrix cl1, SimpleCLMatrix cl2) {
    if (cl1.cellsIdx != cl2.cellsIdx || cl1.lociIdx != cl2.lociIdx)
      throw new RuntimeException
  }
  
  static def void main(String [] args) {
    val matrix = fromCSV(new File(args.get(0)))
    println("nCells: " + matrix.cells.size)
    println("nLoci: " + matrix.loci.size)
  }
}