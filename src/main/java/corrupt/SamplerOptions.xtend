package corrupt

import blang.inits.DesignatedConstructor

class SamplerOptions {
  
  private static SamplerOptions _instance = null
  @DesignatedConstructor
  def static SamplerOptions getInstance() {
    if (_instance === null)
      _instance = new SamplerOptions
    return _instance
  }
  private new() {}
}