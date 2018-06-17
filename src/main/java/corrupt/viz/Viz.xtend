package corrupt.viz

import processing.core.PApplet
import java.io.File

abstract class Viz  {
  
  protected def void draw()
  
  /**
   * Internal coordinate system.
   * width x height
   */
  protected def Pair<? extends Number,? extends Number> size()
  
  protected extension var PApplet applet
  
  var float _declaredWidth = -1
  var float _declaredHeight = -1
  var boolean dimSet = false
  
  /**
   * Width we would like this viz to have, either in pixel as final result, 
   * or when included into another viz. 
   */
  def void declareWidth(Number value) {
    checkNotDeclared
    _declaredWidth = value.floatValue
    _declaredHeight = size.height  * _declaredWidth / size.width
  }
  /**
   * Height we would like this viz to have, either in pixel as final result, 
   * or when included into another viz. 
   */
  def void declareHeight(Number value) {
    checkNotDeclared
    _declaredHeight = value.floatValue
    _declaredWidth = size.width  * _declaredHeight / size.height
  }
  /**
   * When height is set, width is computed automatically from this 
   * Viz internal aspect ratio. Return that inferred width. 
   */
  def float inferWidth() {
    if (!dimSet)
      throw new RuntimeException
    return _declaredWidth
  }
  /**
   * When width is set, height is computed automatically from this 
   * Viz internal aspect ratio. Return that inferred height. 
   */
  def float inferHeight() {
    if (!dimSet)
      throw new RuntimeException
    return _declaredHeight
  }
  
  private def checkNotDeclared() {
    if (dimSet)
      throw new RuntimeException("Can only set one of height, width, the other will be set dynamically")
    dimSet = true
  }
  
  // Use those to avoid confusion in feature calls
  private def void drawViz() { draw }
  private boolean running = false
  
  protected def void addChild(Viz child, Number x, Number y) {
    child.applet = this.applet
    pushMatrix
    pushStyle
    translate(x.floatValue, y.floatValue)
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
    if (!dimSet)
      throw new RuntimeException("Need to set dimension")
    running = true
    val toFile = output !== null
    applet = new PApplet {
      override settings() {
        if (toFile)
          size(_declaredWidth as int, _declaredHeight as int, PDF, output.absolutePath)
        else
          size(_declaredWidth as int, _declaredHeight as int)
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
  
  private def float scaleRatio() { _declaredWidth / size.width }
  
  public static def float width(Pair<? extends Number,? extends Number> p) { return p.key.floatValue }
  public static def float height(Pair<? extends Number,? extends Number> p) { return p.value.floatValue }
}