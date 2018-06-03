package corrupt

import org.junit.Test
import static corrupt.CorruptUtils.syntheticLoci
import blang.engines.internals.factories.Exact
import blang.runtime.SampledModel
import blang.inits.experiments.ExperimentResults
import java.io.File

class UniformTest {
  @Test
  def void testUniformNormalization() {
    val phylo = new PerfectPhylo(1, syntheticLoci(2))
    val model = new Uniform.Builder().setPhylo(phylo).build
    val sModel = new SampledModel(model)
    val exact = new Exact => [
      sampledModel = sModel
      checkLawsGenerateAgreement = true
      results = new ExperimentResults(new File(".")) 
    ]
    exact.performInference
  }
}