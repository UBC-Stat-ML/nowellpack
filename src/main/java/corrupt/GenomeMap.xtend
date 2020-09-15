package corrupt

import java.util.List
import java.util.TreeSet
import java.util.Map
import java.util.LinkedHashMap
import briefj.BriefMaps
import java.util.Comparator
import java.util.ArrayList
import java.util.Collections
import java.util.Collection
import java.util.LinkedHashSet
import org.eclipse.xtend.lib.annotations.Data
import corrupt.post.CLMatrixUtils
import java.io.File

class GenomeMap {
  
  def static List<Locus> orderedLoci(Collection<Locus> _loci) {
    val lociSet = sanitize(_loci)
    if (genomeMapFormatted(lociSet))
      return new GenomeMap(lociSet).orderedLoci
    else
      return lociSet
  }
  
  static def  <T> List<T> sanitize(Collection<T> items) {
    val asSet = new LinkedHashSet(items)
    val asList = new ArrayList(asSet)
    Collections::sort(asList, Comparator.comparing[it.toString])
    return asList  
  }
  
  def static boolean genomeMapFormatted(Collection<Locus> loci) {
    for (locus : loci) 
      if (!genomeMapFormatted(locus))
        return false
    return true
  }
  
  def static boolean genomeMapFormatted(Locus locus) {
    try {
      new ParsedLocus(locus)
    } catch (Exception e) {
      return false
    }
    return true
  }
  
  /**
   * Warning: do not expose the TreeSet as the ordering assumes comparison
   * will only be made with loci in same chromosomes.
   * That's why it's now converted to ArrayList in orderedLoci
   */
  val Map<Integer,TreeSet<Locus>> byChromosome = new LinkedHashMap
  
  new (Collection<Locus> loci) {
    for (locus : loci) {
      val parsed = new ParsedLocus(locus)
      val current = BriefMaps.getOrPut(
        byChromosome, 
        parsed.chr, 
        new TreeSet<Locus>(
          Comparator::comparing[
            val p = new ParsedLocus(it)
            p.leftOneIndexedIncl + p.indexInLocus
          ]
        )
      )
      current.add(locus)
    }
  }
    
  /**
   * When a locus is a difference between consecutive bin, we encode it with 
   * the left locus. As such, binSize is always for a single bin, not spanning 
   * 2 or more loci. 
   */
  @Data
  static class ParsedLocus {
    val int chr
    val String chrString
    val int leftOneIndexedIncl 
    val int rightOneIndexedIncl
    val int indexInLocus
    new (Locus locus) {
      val parsed = locus.toString.split("_")
      if (parsed.size < 4 || parsed.size >5) throw new RuntimeException
      chrString = parsed.get(1)
      chr = chromosomeIndex(chrString)
      leftOneIndexedIncl = Integer.parseInt(parsed.get(2))
      rightOneIndexedIncl = Integer.parseInt(parsed.get(3))
      indexInLocus = if (parsed.size > 4) Integer.parseInt(parsed.get(4)) else 0
    }
    def int binSize() {
      return rightOneIndexedIncl - leftOneIndexedIncl + 1
    }
  }
  
  def static boolean isAdjacent(Locus before, Locus after) {
    val l0 = new ParsedLocus(before)
    val l1 = new ParsedLocus(after)
    if (l0.chr !== l1.chr) return false
    return l0.rightOneIndexedIncl + 1 === l1.leftOneIndexedIncl
  }
  
  def static Locus locus(String chr, int left, int right, int index) {
    return new Locus(chr + "_" + left + "_" + right + (if (index == 0) "" else "_" + index))
  }
  
  def static Locus locus(String chr, int left, int right) {
    locus(chr, left, right, 0)
  } 
  
  def static int chromosomeIndex(String str) {
    val chrStr = str.toUpperCase
    val result = if (chrStr == "X")
      23
    else if (chrStr == "Y")
      24
    else
      Integer.parseInt(chrStr)
    if (result < 1 || result > 24) 
      throw new RuntimeException("Using 1-indexed chromosome.")
    return result
  }
  
  def ArrayList<Locus> neighbors(Locus locus, int neighborhoodSize) {
    val result = new ArrayList<Locus>
    val chrLoci = byChromosome.get(new ParsedLocus(locus).chr)
    var Locus cur = null
    if (neighborhoodSize > 0) {
      cur = locus; for (i : 1..neighborhoodSize) { if (cur !== null) cur = chrLoci.higher(cur); if (cur !== null) result.add(cur) }
      cur = locus; for (i : 1..neighborhoodSize) { if (cur !== null) cur = chrLoci.lower(cur);  if (cur !== null) result.add(cur) }
    }
    return result
  }
  
  def List<Integer> orderedChromosomes() {
    val List<Integer> sortedChrs = new ArrayList(byChromosome.keySet)
    Collections::sort(sortedChrs)
    return sortedChrs
  }
  
  def boolean lociAdjacent() {
    for (chr : orderedChromosomes) {
      val list = new ArrayList(byChromosome.get(chr))
      for (i : 0 ..< (list.size - 1))
        if (!isAdjacent(list.get(i), list.get(i+1)))
          return false
    }
    return true
  }
  
  def List<Locus> orderedLoci(Integer chr) {
    return new ArrayList(byChromosome.get(chr))
  }
  
  def List<Locus> orderedLoci() {
    val List<Locus> result = new ArrayList
    for (chr : orderedChromosomes())
      for (locus : byChromosome.get(chr))
        result.add(locus)
    return result
  }
  
  def static prettyPrintChr(Integer code) {
    if (code == 23) return "X"
    if (code == 24) return "Y"
    if (code < 1 || code > 24)  throw new RuntimeException
    return "" + code
  }
  
  def static void main(String [] args) {
    val matrix = CLMatrixUtils::fromCSV(new File(args.get(0)))
    val map = new GenomeMap(matrix.loci)
    for (chr : map.orderedChromosomes) {
      println(prettyPrintChr(chr) + " " + map.orderedLoci(chr).join(" "))
    }
  }
}