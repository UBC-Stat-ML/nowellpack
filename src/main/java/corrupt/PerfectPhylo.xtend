package corrupt

import java.util.Map
import java.util.LinkedHashMap
import java.util.Set

import static extension corrupt.CorruptExtensionUtils.*
import static corrupt.CorruptStaticUtils.*
import org.eclipse.xtend.lib.annotations.Data
import bayonet.distributions.Random
import org.jgrapht.UndirectedGraph
import briefj.collections.UnorderedPair
import java.util.ArrayList
import java.util.List
import conifer.io.newick.NewickParser
import blang.inits.Input
import blang.inits.DesignatedConstructor
import java.util.LinkedHashSet
import briefj.collections.Tree
import briefj.BriefIO
import java.io.File

@Data class PerfectPhylo {
  val DirectedTree<TreeNode> tree 
  val Set<Locus> loci
  val Set<Cell> cells
  
  def public Map<Cell,Boolean> getTips(Locus locus) { 
      val result = new LinkedHashMap
      _getTips(locus, result, tree.root, false)
      return result
    }
  def private void _getTips(Locus locus, Map<Cell,Boolean> result, TreeNode node, boolean active) {
    switch node {
      Cell : {
        if (!tree.isLeaf(node)) throw new RuntimeException
        result.put(node, active)
      }
      default : {
        val childrenActive = active || node == locus
        for (child : tree.children(node)) 
          _getTips(locus, result, child, childrenActive) 
      }
    }
  }
  
  /**
   * Initialized with a star tree.
   */
  new(Set<Cell> cells, Set<Locus> loci) { 
    this.cells = cells
    this.loci = loci
    tree = new DirectedTree(root)
    for (cell : cells)
      tree.addEdge(root, cell)
    for (locus : loci) 
      tree.addEdge(root, locus)
  }
  
  @DesignatedConstructor
  new(@Input(formatDescription = "Newick string (or 'file ' followed by file from which to load such string)") String newickString) {
    val parser = new NewickParser(
      if (newickString.startsWith("file "))
        BriefIO::fileToString(new File(newickString.replaceFirst("files\\s+", "")))  
      else 
        newickString
    )
    this.cells = new LinkedHashSet
    this.loci = new LinkedHashSet
    this.tree = new DirectedTree(root)
    val Tree<conifer.TreeNode> parseTree = parser.parse
    readParseTree(null, parseTree) 
  }
  
  private def void readParseTree(TreeNode parent, Tree<conifer.TreeNode> parseTree) {
    val parsedNode = parse(parseTree.label.toString)
    if (parent !== null)
      tree.addEdge(parent, parsedNode)
    switch (parsedNode) {
      Cell : cells.add(parsedNode)
      Locus : loci.add(parsedNode)
    }
    for (child : parseTree.children)
      readParseTree(parsedNode, child) 
  }
  
  def void sampleUniform(Random rand) {
    for (cell : cells)
      tree.collapseEdge(cell)
    for (locus : loci)
      tree.collapseEdge(locus)
    val lociTopology = rand.uniformUndirectedTree(loci.size + 1)
    val orderedNodes = new ArrayList<TreeNode>(loci)
    orderedNodes.add(root)
    addSampledLociTopology(lociTopology, loci.size, -1, orderedNodes)
    for (cell : cells)
      tree.addEdge(rand.uniformElement(orderedNodes), cell)
  }
  
  def private void addSampledLociTopology(
    UndirectedGraph<Integer, UnorderedPair<Integer, Integer>> graph, 
    int current, 
    int parent, 
    List<TreeNode> orderedNodes
  ) {
    val currentNode = if (current === loci.size) root else orderedNodes.get(current)
    for (edge : graph.edgesOf(current)) {
      val otherEnd = if (edge.first == current) edge.second else edge.first
      if (otherEnd != parent) {
        val childNode = orderedNodes.get(otherEnd)
        tree.addEdge(currentNode, childNode)
        addSampledLociTopology(graph, otherEnd, current, orderedNodes)
      }
    }
  }
  
  def String toNewick() {
    val result = new StringBuilder()
    toNewick(tree.root, result)
    result.append(";")
    return result.toString
  }
  
  private def void toNewick(TreeNode node, StringBuilder builder) {
    val children = tree.children(node)
    if (!children.empty) {
      builder.append("(")
      for (var int cIndex = 0; cIndex < children.size(); cIndex++) {
        toNewick(children.get(cIndex), builder)
        if (cIndex !== children.size - 1)
          builder.append(",")
      }
      builder.append(")")
    }
    val label = node.toString
    if (label.contains("(") || 
        label.contains(")") || 
        label.contains(",") || 
        label.contains(":") ||
        label.contains(";"))
      throw new RuntimeException();
    if (node != root) 
      builder.append(label) 
  }
  
  override String toString() { toNewick }
}