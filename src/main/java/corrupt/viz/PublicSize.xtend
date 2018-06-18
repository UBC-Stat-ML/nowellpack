package corrupt.viz

import org.eclipse.xtend.lib.annotations.Data
import blang.inits.DesignatedConstructor
import blang.inits.Input
import java.util.List

/**
 * When a viz is used by something else (either a main application 
 * that will print the viz to screen/pdf, or another viz), either a
 * height or width needs to be specified, but not both since the ratio 
 * is already determined by the PrivateSize.
 */
@Data
class PublicSize {
  val boolean isHeight
  val float value
  
  public static val format = "(height|width) <value>"
  @DesignatedConstructor
  static def PublicSize parseSize(@Input(formatDescription = format) List<String> arguments) {
    if (arguments.size !== 2)
      throw new RuntimeException(arguments.join(" ") + " does not comply " + format)
    val isHeight = 
      switch (arguments.get(0).toLowerCase.trim) {
        case "height" : true
        case "width" : false
        default : throw new RuntimeException(arguments.join(" ") + " does not comply " + format)
      }
    val value = Float.parseFloat(arguments.get(1))
    return new PublicSize(isHeight, value)
  }
}

  