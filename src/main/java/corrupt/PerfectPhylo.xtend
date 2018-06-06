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
  new(@Input String newickString) {
    val parser = new NewickParser(newickString)
    this.cells = new LinkedHashSet
    this.splits = new LinkedHashMap
    this.tree = new DirectedTree(root)
    val Tree<conifer.TreeNode> parseTree = parser.parse
    setTreeFrom(parseTree) 
  }
  
  private def setTreeFrom(Tree<conifer.TreeNode> parseTree) {
    val parsedNode = parse(parseTree.label.toString)
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
    toNewick(tree.root, result, new ArrayList)
    result.append(";")
    return result.toString
  }
  
  public static val COLLAPSED_LOCI_SEP = "+"
  private def void toNewick(TreeNode node, StringBuilder builder, List<String> loci) {
    val children = tree.children(node)
    if (node != root)
        loci.add(node.toString)
    if (children.size == 1 && children.get(0) instanceof Locus) {
      // collapse lists of loci
      toNewick(children.get(0), builder, loci)
    } else {
      if (!children.empty) {
        builder.append("(")
        for (var int cIndex = 0; cIndex < children.size(); cIndex++) {
          toNewick(children.get(cIndex), builder, new ArrayList)
          if (cIndex !== children.size - 1)
            builder.append(",")
        }
        builder.append(")")
      }
      val label = loci.join(COLLAPSED_LOCI_SEP)
      if (label.contains("(") || 
          label.contains(")") || 
          label.contains(",") || 
          label.contains(":") ||
          label.contains(";"))
        throw new RuntimeException();
      builder.append(label);
    }
  }
  
  override String toString() { toNewick }
}