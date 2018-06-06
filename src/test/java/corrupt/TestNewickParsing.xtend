package corrupt

import org.junit.Test



import static extension corrupt.CorruptExtensionUtils.*
import static corrupt.CorruptStaticUtils.*
import static org.junit.Assert.*
import bayonet.distributions.Random

class TestNewickParsing {
  @Test
  def testRoundTrip() {
    val tree = new PerfectPhylo(syntheticCells(10), syntheticLoci(10))
    tree.sampleUniform(new Random(1))
    val newick = tree.toNewick
    println(newick)
    val tree2 = new PerfectPhylo(newick)
    assertEquals(tree, tree2)
  }
}
