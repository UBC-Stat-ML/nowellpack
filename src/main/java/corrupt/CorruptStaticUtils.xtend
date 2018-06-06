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
  
  def static List<TreeNode> parse(String description) {
         if (description.startsWith(Locus::PREFIX)) return                                 loci(description.replaceFirst(Locus::PREFIX, ""))
    else if (description.startsWith(Cell::PREFIX )) return Collections::singletonList(new Cell (description.replaceFirst(Cell::PREFIX,  "")))
    else if (description == root.toString) return Collections::singletonList(root)
    else throw new RuntimeException
  }
  
  private def static List<TreeNode> loci(String string) {
    throw new UnsupportedOperationException("TODO: auto-generated method stub")
  }
}
