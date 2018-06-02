package corrupt

import corrupt.DirectedTree
import org.junit.Test
import static org.junit.Assert.assertTrue

class TestTree {
  @Test
  def void testInverse() {
    val t = smallTree
    val movedChildren = t.collapseEdge(1)
    assertTrue(t != smallTree)
    t.addEdge(0, 1, movedChildren)
    assertTrue(t == smallTree)
  }
  
  def DirectedTree<Integer> smallTree() {
    val result = new DirectedTree(0) => [
      addEdge(0, 1)
        addEdge(1, 5)
        addEdge(1, 6)
      addEdge(0, 2)
      addEdge(0, 3)
      addEdge(0, 4)
    ]
    return result
  }
}