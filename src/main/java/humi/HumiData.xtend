package humi

import blang.io.GlobalDataSource
import blang.inits.Arg
import blang.types.Plate
import blang.types.Plated
import com.rits.cloning.Immutable
import blang.types.Index
import blang.inits.DefaultValue

@Immutable
class HumiData {
  @Arg public GlobalDataSource source
  
  @Arg public Plated<CountFrequencies> histograms
  
  @Arg public Plate<Integer> targets
  @Arg public Plate<String> genes
  @Arg public Plate<String> experiments
  
  @Arg            @DefaultValue("ctl")
  public String controlPrefix = "ctl"
  def boolean isControl(Index<String> gene) {
    if (gene.plate != genes) throw new RuntimeException
    return gene.key.startsWith(controlPrefix)
  }
}