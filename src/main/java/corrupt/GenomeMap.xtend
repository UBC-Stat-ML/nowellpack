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

class GenomeMap {
  
  def static List<Locus> orderLoci(Collection<Locus> _loci) {
    val lociSet = sanitize(_loci)
    if (genomeMapFormatted(lociSet))
      return new GenomeMap(lociSet).orderLoci
    else
      return lociSet
  }
   
  def static List<Locus> orderLoci(Collection<Locus> _loci, String pt) {
    val lociSet = sanitize(_loci)
    val lociSet_copy = sanitize(_loci)
    for (locus : lociSet_copy){
    	  if (locus.getPrintType() != pt){
    	  	lociSet.remove(locus)
    	  }
    }
    if (genomeMapFormatted(lociSet))
      return new GenomeMap(lociSet).orderLoci
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
  
  val Map<Integer,TreeSet<Locus>> byChromosome = new LinkedHashMap
  
  new (Collection<Locus> loci) {
    for (locus : loci) {
      val parsed = new ParsedLocus(locus)
      val current = BriefMaps.getOrPut(byChromosome, parsed.chr, new TreeSet<Locus>(Comparator::comparing[Locus l | new ParsedLocus(l).leftOneIndexedIncl]))
      current.add(locus)
    }
  }
  
  private static class ParsedLocus {
    val int chr
    val int leftOneIndexedIncl
    val int rightOneIndexedIncl
    new (Locus locus) {
      val parsed = locus.toString.split("_")
      val chrStr = parsed.get(1).toUpperCase
      if (chrStr == "X")
        chr = 23
      else if (chrStr == "Y")
        chr = 24
      else
        chr = Integer.parseInt(chrStr)
      leftOneIndexedIncl = Integer.parseInt(parsed.get(2))
      rightOneIndexedIncl = Integer.parseInt(parsed.get(3))
    }
  }
  
  def ArrayList<Locus> neighbors(Locus locus, int neighborhoodSize) {
    val result = new ArrayList<Locus>
    val chrLoci = byChromosome.get(new ParsedLocus(locus).chr)
    var Locus cur = null
    cur = locus; for (i : 1..neighborhoodSize) { if (cur !== null) cur = chrLoci.higher(cur); if (cur !== null) result.add(cur) }
    cur = locus; for (i : 1..neighborhoodSize) { if (cur !== null) cur = chrLoci.lower(cur);  if (cur !== null) result.add(cur) }
    return result
  }
  
  def List<Integer> orderedChromosomes() {
    val List<Integer> sortedChrs = new ArrayList(byChromosome.keySet)
    Collections::sort(sortedChrs)
    return sortedChrs
  }
  
  def Collection<Locus> orderedLoci(Integer chr) {
    return byChromosome.get(chr)
  }
  
  def List<Locus> orderLoci() {
    val List<Locus> result = new ArrayList
    for (chr : orderedChromosomes())
      for (locus : byChromosome.get(chr))
        result.add(locus)
    return result
  }
}