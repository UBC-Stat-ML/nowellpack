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
    .setOverdispersion(StaticUtils::latentReal)
      .build 
      
  val MarginalizedSingleGuide margSingleGuideModel = new MarginalizedSingleGuide.Builder()
    .setFrequencies(new SimpleCountFrequencies)
    .setMean(StaticUtils::latentReal) 
    .setShape(StaticUtils::latentReal)
    .setRate(StaticUtils::latentReal)
    .setOverdispersion(StaticUtils::latentReal)
      .build 
  
  @Test 
  def void exactInvarianceTest() {
    val test = new ExactInvarianceTest
    test.nPosteriorSamplesPerIndep = 500 
    test.add(new Instance(singleGuideModel, 
      [mean.doubleValue], 
      [overdispersion.doubleValue], 
      [NUMIMean.doubleValue]
    ))  
    test.add(new Instance(margSingleGuideModel, 
      [mean.doubleValue], 
      [overdispersion.doubleValue], 
      [shape.doubleValue],
      [rate.doubleValue]
    ))  
    test.check
  }
  
}