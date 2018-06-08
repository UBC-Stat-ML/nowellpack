package corrupt

import java.util.Set
import java.util.List
import java.util.Collections

class CorruptStaticUtils {
  public static val TreeNode root = new TreeNode("ROOT")
  
  def static Set<Locus> syntheticLoci(int n) {
    (0 ..< n).map[new Locus(it.toString)].toSet
  }
  
  def static Set<Cell> syntheticCells(int n) {
    (0 ..< n).map[new Cell(it.toString)].toSet
  }
  
  def static double logNPerfectPhylo(int nCells, int nLoci) {
    return (nLoci + nCells - 1) * Math.log(nLoci + 1)
  }
  
  def static TreeNode parse(String description) {
    // Note: taking care of removing the prefix is taken care of in constructors
         if (description.startsWith(Locus::PREFIX)) return new Locus(description)
    else if (description.startsWith(Cell::PREFIX )) return new Cell (description)
    else return root
  }
}