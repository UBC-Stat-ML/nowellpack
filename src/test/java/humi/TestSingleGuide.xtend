package humi

import org.junit.Test
import blang.validation.ExactInvarianceTest
import blang.validation.Instance
import blang.types.StaticUtils

class TestSingleGuide {
  
  val SingleGuide singleGuideModel = new SingleGuide.Builder()
    .setFrequencies(new SimpleCountFrequencies)
    .setMean(StaticUtils::latentReal)
    .setNUMIMean(StaticUtils::latentReal)
//    .setOverdispersion(StaticUtils::latentReal)
      .build 
  
  @Test 
  def void exactInvarianceTest() {
    val test = new ExactInvarianceTest
    test.nPosteriorSamplesPerIndep = 500
    test.add(new Instance(singleGuideModel, 
      [mean.doubleValue], 
//      [overdispersion.doubleValue], 
      [NUMIMean.doubleValue]
    ))  
    test.check
  }
  
}