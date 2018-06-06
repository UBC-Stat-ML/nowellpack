package corrupt.viz

import processing.core.PApplet
import java.io.File

abstract class Viz  {
  
  protected def void draw()
  
  /**
   * width x height
   */
  protected def Pair<Float,Float> size()
  
  protected extension var PApplet applet
  
  public var float width = 500 // adjust height auto
  
  // Use those to avoid confusion in feature calls
  private def float declaredWidth() { this.width }
  private def void drawViz() { draw }
  private boolean running = false
  
  protected def void addChild(Viz child, float x, float y) {
    child.applet = this.applet
    pushMatrix
    pushStyle
    translate(x, y)
    scale(child.scaleRatio, child.scaleRatio)
    child.drawViz
    popMatrix
    popStyle
  }
  
  def void show() { execute(null) }
  def void output(String path) { output(new File(path)) }
  def void output(File output) { execute(output) }
  private def void execute(File output) {
    if (running)
      throw new RuntimeException("Cannot run twice.")
    running = true
    val toFile = output !== null
    val size = size()
    val canvaHeight = size.height  * this.width / size.width
    applet = new PApplet {
      override settings() {
        if (toFile)
          size(declaredWidth() as int, canvaHeight as int, PDF, output.absolutePath)
        else
          size(declaredWidth() as int, canvaHeight as int)
      }
      override draw() {
        scale(scaleRatio, scaleRatio)
        drawViz()
        if (toFile)
          exit
      }
      override exitActual() {
        if (!toFile) 
          super.exitActual
      }
    }
    PApplet.runSketch(#[this.class.simpleName.toString], applet) 
  }
  
  private def float scaleRatio() { declaredWidth() / size.width }
  
  private def float width(Pair<Float,Float> p) { return p.key }
  private def float height(Pair<Float,Float> p) { return p.value }
}