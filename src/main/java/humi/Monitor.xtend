package humi

import blang.core.RealVar

class Monitor implements RealVar {
  var RealVar monitored = null
  
  new(RealVar monitored) {
    this.monitored = monitored
  }
  
  // need this to use as plated
  // then need to use init
  new() { 
    monitored = null
  }
  
  def void init(RealVar monitored) {
    if (this.monitored !== null) 
      throw new RuntimeException("Already monitored")
    this.monitored = monitored
  }
  
  override doubleValue() {
    if (monitored === null) throw new RuntimeException
    return monitored.doubleValue
  }
  
  override toString() { return "" + doubleValue }
}