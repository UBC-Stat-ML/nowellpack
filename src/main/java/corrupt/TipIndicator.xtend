package corrupt

import org.eclipse.xtend.lib.annotations.Accessors

class TipIndicator {
  @Accessors(PUBLIC_GETTER, PACKAGE_SETTER)
  var boolean included = false // Changes here only triggered by topology changes. Do not change manually.
  @Accessors(PUBLIC_GETTER) val Cell cell
  
  new (Cell cell) {
    this.cell = cell
  }
  
  override toString() { "tip(" + included + ")" }
}