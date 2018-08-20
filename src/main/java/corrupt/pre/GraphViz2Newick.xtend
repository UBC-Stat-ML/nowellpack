package corrupt.pre

import blang.inits.Arg
import java.io.File
import blang.inits.experiments.Experiment
import java.util.ArrayList
import corrupt.TreeNode
import briefj.BriefIO
import briefj.BriefStrings
import java.util.regex.Pattern
import corrupt.CorruptStaticUtils
import corrupt.Locus
import corrupt.Cell
import corrupt.PerfectPhylo
import java.util.LinkedHashSet
import corrupt.DirectedTree
import org.jgrapht.DirectedGraph
import bayonet.graphs.GraphUtils

class GraphViz2Newick extends Experiment {
  @Arg 
  public File input
  
  val edgePattern = Pattern::compile("([^-]+)[ ][-][>][ ]([^-]+);")
  
  override run() {
    val edges = new ArrayList<Pair<TreeNode, TreeNode>>
    for (line : BriefIO::readLines(input)) 
      if (line.contains("->")) {
        val split = BriefStrings::allGroupsFromFirstMatch(edgePattern, line)
        edges.add(parse(split.get(0)) -> parse(split.get(1)))
    }
    val cells = new LinkedHashSet<Cell> 
    val loci = new LinkedHashSet<Locus>
    val iterator = edges.iterator
    while (iterator.hasNext) {
      val edge = iterator.next
      for (endPt : #[edge.key, edge.value]) 
        switch (endPt) {
          Locus : loci.add(endPt)
          Cell : {
            // SCITE puts cells in duplicates sometimes, remove these duplicates
            if (cells.contains(endPt))
              iterator.remove
            else
              cells.add(endPt)
          }
        }
      }
    val DirectedGraph<TreeNode, org.apache.commons.lang3.tuple.Pair<TreeNode,TreeNode>> topo = GraphUtils.<TreeNode>newDirectedGraph
    for (edge : edges) {
      topo.addVertex(edge.key)
      topo.addVertex(edge.value)
      topo.addEdge(edge.key, edge.value)
    }
    val tree = new DirectedTree(CorruptStaticUtils::root, topo)
    val perfectPhylo = new PerfectPhylo(cells, loci, tree)
    val outputStr = perfectPhylo.toNewick
    BriefIO::write(results.getFileInResultFolder("output.newick"), outputStr) 
  }
  
  def TreeNode parse(String str) {
    if (str == "Root") return CorruptStaticUtils::root
    else if (str.startsWith(Locus::PREFIX)) return new Locus(str)
    else return new Cell(str.substring(1)) 
  }
  
  def public static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
}