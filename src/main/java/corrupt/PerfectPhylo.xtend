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
    for (locus : loci) 
      splits.put(locus, Split::initializeEmpty(tree, locus, cells))
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
  
  override String toString() { tree.toString }
}