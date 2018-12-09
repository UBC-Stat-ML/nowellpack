package corrupt

import blang.inits.Arg
import blang.inits.DefaultValue
import blang.inits.DesignatedConstructor

class SamplerOptions {
  @Arg                      @DefaultValue("true")
  public boolean useCellReallocationMove = true
  
  private static SamplerOptions _instance = null
  @DesignatedConstructor
  def static SamplerOptions getInstance() {
    if (_instance === null)
      _instance = new SamplerOptions
    return _instance
  }
  private new() {}
}