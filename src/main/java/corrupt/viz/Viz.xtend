package corrupt.viz

import processing.core.PApplet
import java.io.File

abstract class Viz  {
  def void draw(extension PApplet applet)
  
  /**
   * width x height
   */
  def Pair<Float,Float> size()
  
  public var float canvaWidth = 500 // adjust height auto
  
  def void show() { execute(null) }
  def void output(String path) { output(new File(path)) }
  def void output(File output) { execute(output) }
  private def void execute(File output) {
    val toFile = output !== null
    val size = size()
    val canvaHeight = size.height  * canvaWidth / size.width
    val applet = new PApplet {
      override settings() {
        if (toFile)
          size(canvaWidth as int, canvaHeight as int, PDF, output.absolutePath)
        else
          size(canvaWidth as int, canvaHeight as int)
      }
      override draw() {
        val scaleRatio = canvaWidth / size.width 
        scale(scaleRatio, scaleRatio)
        draw(this)
        if (toFile)
          exit
      }
      override exitActual() {}
    }
    PApplet.runSketch(#[this.class.simpleName.toString], applet) 
  }
  
  private def float width(Pair<Float,Float> p) { return p.key }
  private def float height(Pair<Float,Float> p) { return p.value }
}