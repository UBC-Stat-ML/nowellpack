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
import java.util.Collections

@Data class PerfectPhylo {
  val DirectedTree<TreeNode> tree 
  val Set<Locus> loci
  val Set<Cell> cells 
  
  def Map<Cell,Boolean> getTips(Locus locus) { 
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
    this(cells, loci, new DirectedTree(root))
    for (cell : cells)
      tree.addEdge(root, cell)
    for (locus : loci) 
      tree.addEdge(root, locus)
  }
  
  new(Set<Cell> cells, Set<Locus> loci, DirectedTree<TreeNode> tree) {
    this.cells = cells
    this.loci = loci
    this.tree = tree
  }
  
  def static PerfectPhylo generateUniform(int nCells, int nLoci, Random rand) {
    val result = new PerfectPhylo(CorruptStaticUtils::syntheticCells(nCells), CorruptStaticUtils::syntheticLoci(nLoci))
    result.sampleUniform(rand)
    return result
  }
  
  def static PerfectPhylo parseNewick(File f) {
    return new PerfectPhylo(BriefIO.fileToString(f)) 
  }
  
  @DesignatedConstructor
  new(@Input(formatDescription = "Newick string (or 'file ' followed by file from which to load such string)") String newickString) {
    val parser = new NewickParser(
      if (newickString.startsWith("file "))
        BriefIO::fileToString(new File(newickString.replaceFirst("file\\s+", "")))  
      else 
        newickString
    )
    this.cells = new LinkedHashSet
    this.loci = new LinkedHashSet
    this.tree = new DirectedTree(root)
    val Tree<conifer.TreeNode> parseTree = parser.parse
    if (parse(parseTree.label.toString) != root)
      throw new RuntimeException("PerfectPhylo newick file should have root called " + root)
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
  
  def void set(PerfectPhylo phylo) {
    if (this.cells != phylo.cells || this.loci != phylo.loci) throw new RuntimeException
    for (cell : cells)
      tree.collapseEdge(cell)
    for (locus : loci)
      tree.collapseEdge(locus)
    for (edge : phylo.tree.edges) 
      this.tree.addEdge(edge.left, edge.right)
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
  
  /**
   * A tree constructed from the perfect phylo by: 
   * 1) removing all subtrees that contain no cell, and,
   * 2) collapse chains into a single set containing all the nodes in the chain (which may then be a mix of cells, loci, root).
   */
  def DirectedTree<Set<TreeNode>> collapsedTree() {
    val collapsed = collapsedTree(tree.root)
    val result = new DirectedTree(collapsed.label)
    copy(collapsed, result)
    return result
  }
  
  private static def void copy(Tree<Set<TreeNode>> src, DirectedTree<Set<TreeNode>> dest) {
    for (child : src.children) {
      dest.addEdge(src.label, child.label)
      copy(child, dest)
    }
  }
  
  private def Tree<Set<TreeNode>> collapsedTree(TreeNode node) {
    val recursions = new ArrayList<Tree<Set<TreeNode>>>
    for (child : tree.children(node)) {
      val current = collapsedTree(child)
      if (current !== null)
        recursions.add(current)
    }
    val nodeSet = new LinkedHashSet(Collections.singleton(node))
    if (recursions.size === 0) {
      if (node instanceof Cell)
        return new Tree(nodeSet)
      else
        return null
    } else if (recursions.size === 1) {
      val child = recursions.get(0)
      if (child === null)
        throw new RuntimeException
      child.label.add(node) 
      return child
    } else {
      return new Tree(nodeSet, recursions)
    }
  }
  
  override String toString() { toNewick }
}