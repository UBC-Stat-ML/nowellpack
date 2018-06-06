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
  // Warning: updates on the tree need to be mirrored to the splits
  val DirectedTree<TreeNode> tree 
  val Map<Locus, Split> splits
  val Set<Cell> cells
  
  /**
   * Initialized with a star tree.
   */
  new(Set<Cell> cells, Set<Locus> loci) { 
    this.cells = cells
    tree = new DirectedTree(root)
    splits = new LinkedHashMap
    for (cell : cells)
      tree.addEdge(root, cell)
    for (locus : loci) {
      splits.put(locus, Split::initializeEmpty(tree, locus, cells))
      tree.addEdge(root, locus)
    }
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
    this.splits = new LinkedHashMap
    this.tree = new DirectedTree(root)
    val Tree<conifer.TreeNode> parseTree = parser.parse
    readParseTree(null, parseTree) 
    for (locus : loci)
      splits.put(locus, Split::initializeEmpty(tree, locus, cells))
    updateAllSplits
  }
  
  private def void readParseTree(TreeNode parent, Tree<conifer.TreeNode> parseTree) {
    val parsedNode = parse(parseTree.label.toString)
    if (parent !== null)
      tree.addEdge(parent, parsedNode)
    switch (parsedNode) {
      Cell : cells.add(parsedNode)
      Locus : splits.put(parsedNode, null)
    }
    for (child : parseTree.children)
      readParseTree(parsedNode, child) 
  }
  
  def Set<Locus> getLoci() { splits.keySet }
  
  def TipIndicator tipIndicator(Cell cell, Locus locus) {
    return splits.get(locus).tipIndicators.get(cell)
  }
  
  def void updateAllSplits() {
    for (split : splits.values)
      split.updateTips
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
    updateAllSplits
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
  
  override hashCode() {
    return tree.hashCode
  }
  
  override boolean equals(Object obj) {
    if (this === obj)
      return true
    if (obj === null)
      return false
    if (getClass() !== obj.getClass())
      return false
    val other = obj as PerfectPhylo
    if (this.tree === null) {
      if (other.tree !== null)
        return false
    } else if (!this.tree.equals(other.tree))
      return false
    return true
  }
}