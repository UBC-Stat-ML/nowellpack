package humi

import blang.core.RealVar

// Not ideal yet / ready to move to SDK
// Problem is that lambdas in xtext are not always treated 
// the same: in eclipse, transpiled into anoymous classes, so all good,
// in gradle, instead they become java 8 lambdas, and hence crash deep copying
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
    if (initialized) 
      throw new RuntimeException("Already monitored")
    this.monitored = monitored
  }
  
  def boolean initialized() {
    return this.monitored !== null
  }
  
  override doubleValue() {
    if (monitored === null) 
      throw new RuntimeException
    return monitored.doubleValue
  }
  
  override toString() { return "" + doubleValue }
}