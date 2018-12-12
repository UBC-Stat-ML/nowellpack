package corrupt

import org.junit.Test

import static corrupt.CorruptStaticUtils.*

import bayonet.distributions.Random
import corrupt.post.CLMatrixUtils
import corrupt.post.SimpleCLMatrix
import java.util.LinkedHashSet
import org.junit.Assert
import org.junit.Rule
import org.junit.rules.TemporaryFolder
import blang.inits.experiments.ExperimentResults

class GrowTreeTest {
  
  @Rule
  public TemporaryFolder folder = new TemporaryFolder();
  
  val nCells = 1000
  val nLoci = 1000
  
  val rand = new Random(1)
  
  @Test
  def void testCellPlacement() {
    val cells = syntheticCells(nCells)
    val loci = syntheticLoci(nLoci)
    val phylo = new PerfectPhylo(cells, loci)
    phylo.sampleUniform(rand)
    val tipInclPrs = CLMatrixUtils::syntheticInclusionPrs(rand, phylo, 0.005) 
    
    // remove one cell
    val removed = cells.iterator.next
    val parent = phylo.tree.parent(removed)
    phylo.tree.collapseEdge(removed)
    
    // split the matrix
    val c0 = new LinkedHashSet(cells) => [ remove(removed) ]
    val c1 = new LinkedHashSet() => [add(removed) ]
    val SimpleCLMatrix m0 = new SimpleCLMatrix(c0, loci)
    val SimpleCLMatrix m1 = new SimpleCLMatrix(c1, loci)
    
    for (locus : loci) 
      for (cell : cells) 
        if (cell == removed) {
          m1.set(cell, locus, tipInclPrs.get(cell, locus))
        } else {
          m0.set(cell, locus, tipInclPrs.get(cell, locus))
        }
      
    val r = new GrowTree => [
      results = new ExperimentResults(folder.newFolder)
    ]
    r.matrix = m1
    r.phylo = new PerfectPhylo(c0, loci, phylo.tree)
    r.run    
    
    val reconstructedParent = phylo.tree.parent(removed)
    
    Assert.assertEquals(println(parent), reconstructedParent)
  }
  
  @Test
  def void testLocusPlacement() {
    val cells = syntheticCells(nCells)
    val loci = syntheticLoci(nLoci)
    val phylo = new PerfectPhylo(cells, loci)
    phylo.sampleUniform(rand)
    val tipInclPrs = CLMatrixUtils::syntheticInclusionPrs(rand, phylo, 0.005) 
    
    // remove one cell
    val removed = loci.iterator.next
    val parent = phylo.tree.parent(removed)
    phylo.tree.collapseEdge(removed)
    
    // split the matrix
    val l0 = new LinkedHashSet(loci) => [ remove(removed) ]
    val l1 = new LinkedHashSet() => [add(removed) ]
    val SimpleCLMatrix m0 = new SimpleCLMatrix(cells, l0)
    val SimpleCLMatrix m1 = new SimpleCLMatrix(cells, l1)
    
    for (locus : loci) 
      for (cell : cells) 
        if (locus == removed) {
          m1.set(cell, locus, tipInclPrs.get(cell, locus))
        } else {
          m0.set(cell, locus, tipInclPrs.get(cell, locus))
        }
     
    val r = new GrowTree => [
      results = new ExperimentResults(folder.newFolder)
    ]
    r.matrix = m1
    r.phylo = new PerfectPhylo(cells, l0, phylo.tree)
    r.run    
    
    val reconstructedParent = phylo.tree.parent(removed)
    
    Assert.assertEquals(println(parent), reconstructedParent)
  }
}