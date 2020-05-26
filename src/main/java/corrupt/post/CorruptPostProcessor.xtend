package corrupt.post

import blang.runtime.internals.DefaultPostProcessor
import blang.runtime.Runner
import java.io.File
import blang.inits.parsing.ConfigFile
import corrupt.Greedy
import briefj.BriefIO
import corrupt.PerfectPhylo
import xlinear.Matrix
import static xlinear.MatrixOperations.*
import static extension xlinear.MatrixExtensions.*
import corrupt.CorruptPhylo
import blang.inits.experiments.tabwriters.TabularWriter
import blang.inits.experiments.Experiment
import corrupt.Locus
import java.util.Map
import java.util.LinkedHashMap
import corrupt.viz.PerfectPhyloViz
import viz.core.Viz
import blang.inits.experiments.tabwriters.factories.CSV

import static blang.inits.experiments.tabwriters.factories.CSV.csvFile 

class CorruptPostProcessor extends DefaultPostProcessor  {
  
  File binaryObservationsFile
  
  var File sampleDir
  override run() {
    // load run
    sampleDir = new File(blangExecutionDirectory.get, Runner::SAMPLES_FOLDER)
    val args = ConfigFile::parse(new File(blangExecutionDirectory.get, DETAILED_ARGUMENT_FILE))
    val binaryObsPath = args.child("model").child("binaryMatrix").argumentValue.get.join
    binaryObservationsFile = new File(binaryObsPath)
    val binaryObservations = BinaryCLMatrix::create(binaryObservationsFile.absolutePath)
    
    // postprocessing steps:
    val consensus = decode()
    predictiveResults()
    goodnessOfFit(consensus)
    
    // call viz
    if (runPxviz) 
      treeViz(consensus.reconstruction, binaryObservations)
    
    super.run
  }
  
  def getCsvFile(File f, String s) { CSV::csvFile(f, s) } // needed to sidestep obscure scope corner case
  
  def treeViz(PerfectPhylo reconstruction, BinaryCLMatrix observations) {
    // for the support values, don't reuse the ones from consensus construction; they might use the logistic transform!
    val child = results.child("averagingNoLogisticsTransform")
    val averager = new AverageCLMatrices => [
      csvFile = getCsvFile(sampleDir, "phylo")   
      results = child
      logisticTransform = false
    ]
    averager.run
    val supportValues = ReadOnlyCLMatrix::readOnly(averager.average)
    val convertedObs = ReadOnlyCLMatrix::readOnly(observations)
    val publicSize = Viz.fixWidth(500)
    val treeViz = new PerfectPhyloViz(reconstruction, #[supportValues,convertedObs], publicSize)
    val treeVizResults = results.child("treeViz") 
    treeViz.output(treeVizResults.getFileInResultFolder("consensus.pdf"))
  }
  
  def goodnessOfFit(CorruptPhylo consensus) {
    val dataMatrix = BinaryCLMatrix::create(binaryObservationsFile.absolutePath)
    val gofResults = results.child("goodnessOfFit")
    
    // for the consensus
    val consensusGof = new GoodnessOfFit(dataMatrix, consensus.reconstruction)
    val consensusGofWriter = gofResults.getTabularWriter("consensus")
    consensusGof.logTo(consensusGofWriter)
    
    // 2D plot for the consensus
    val consensusGofWriterByLocus = gofResults.getTabularWriter("consensus-by-locus")
    consensusGof.logToByLocus(consensusGofWriterByLocus)
    
    // for the trace
    val phyloSamples = csvFile(sampleDir, "phylo")
    val traceGofWriter = gofResults.getTabularWriter("trace")
    for (line : BriefIO::readLines(phyloSamples).indexCSV) {
      val phyloStr = line.get("value")
      val tree = new PerfectPhylo(phyloStr)
      val sampleGof = new GoodnessOfFit(dataMatrix, tree)
      sampleGof.logTo(traceGofWriter.child("sample", line.get("sample"))) 
    }
    
    gofResults.closeAll
    
    val script = '''
    require("ggplot2")
    require("dplyr")
    
    trace <- read.csv("«csvFile(gofResults.resultsFolder, traceGofWriter.name)»")
    
    consensus <- read.csv("«csvFile(gofResults.resultsFolder, consensusGofWriter.name)»")
    
    «FOR stat : #["gof", "empiricalFN", "empiricalFP"]»
    p <- ggplot(trace, aes(x = sample, y = «stat»)) + geom_line() + geom_hline(yintercept=consensus$«stat») + theme_bw()
    ggsave("«gofResults.getFileInResultFolder(stat + ".pdf")»", p)
    «ENDFOR»
    
    consensusByLocus <- read.csv("«csvFile(gofResults.resultsFolder, consensusGofWriterByLocus.name)»")
    p <- ggplot(consensusByLocus, aes(x = empiricalFN, y = empiricalFP)) + geom_point() + theme_bw()
    ggsave("«gofResults.getFileInResultFolder("consensus-by-locus.pdf")»", p)
    '''
    callR(gofResults.getFileInResultFolder(".script.r"), script)
  }
  
  def CorruptPhylo decode() {
    val child = results.child("averaging")
    val averager = new AverageCLMatrices => [
      csvFile = getCsvFile(sampleDir, "phylo")
      results = child
    ]
    averager.run
    
    val greedyDecoder = new Greedy => [
      tipInclusionProbabilities = ReadOnlyCLMatrix::create(child.getFileInResultFolder(AverageCLMatrices::OUTPUT_NAME).absolutePath) 
    ]
    val consensus = greedyDecoder.infer
    BriefIO::write(results.getFileInResultFolder("consensus.newick"), consensus.toString)
    return consensus
  }
  
  static class GoodnessOfFit {
    public val Matrix counts = dense(2,2)
    public val Map<Locus,Matrix> byLocus = new LinkedHashMap
    new () {}
    new (BinaryCLMatrix observations, PerfectPhylo reconstruction) {
      for (locus : observations.loci) {
        val locusSpecific = dense(2,2)
        byLocus.put(locus, locusSpecific)
        val tips = reconstruction.getTips(locus)
        for (cell : tips.keySet) {
          for (m : #[counts, locusSpecific])
            m.increment(if (tips.get(cell)) 1 else 0, observations.get(cell, locus) as int, 1.0)
        }
      }
    }
    def static double gof(Matrix counts) { return (counts.get(0,0) + counts.get(1,1)) / counts.sum }
    def static double empiricalFN(Matrix counts) { return counts.get(1, 0) / counts.row(1).sum }
    def static double empiricalFP(Matrix counts) { return counts.get(0, 1) / counts.row(0).sum }
    def logTo(TabularWriter writer) {
      writer.write(
        "gof" -> gof(counts),
        "empiricalFN" -> empiricalFN(counts),
        "empiricalFP" -> empiricalFP(counts)
      )
    }
    def logToByLocus(TabularWriter writer) {
      for (locus : byLocus.keySet)
        writer.write(
          "locus" -> locus,
          "gof" -> gof(byLocus.get(locus)),
          "empiricalFN" -> empiricalFN(byLocus.get(locus)),
          "empiricalFP" -> empiricalFP(byLocus.get(locus))
        )
    }
  }
  
  def predictiveResults() {
    val outputFolder = results.child("predictivePlots")
    val predictiveTraces = csvFile(sampleDir, "predictives") 
    if (BriefIO::readLines(predictiveTraces).head === null) return
    
    val script = '''
    require("ggplot2")
    require("dplyr")
    library("reshape2")
    
    # load data
    
    dataRaw <- read.csv("«binaryObservationsFile»") 
        names(dataRaw)[names(dataRaw) == 'cells'] <- 'cell'
        names(dataRaw)[names(dataRaw) == 'loci'] <- 'locus'
    
    # load raw MCMC trace
    
    raw <- read.csv("«predictiveTraces»")
    
    names(raw)[names(raw) == 'map_key_0'] <- 'locus'
    names(raw)[names(raw) == 'map_key_2'] <- 'cell'
    
    # Remove burn in
    
    n_samples <- max(raw$«Runner.sampleColumn»)
    cut_off <- n_samples * «burnInFraction»
    raw <- subset(raw, «Runner.sampleColumn» > cut_off)
    
    
    # Calibration viz   
    
    presencePredictivesTrace <- raw %>% filter(predicted == 'presence')
    
    presencePosteriors <- presencePredictivesTrace %>%
          group_by(locus, cell) %>%
          mutate(weight = 1.0/value) %>%
          summarise(
            simplePresence = mean(value),
            leaveOneOutPresence = 1.0 / mean(weight), # harmonic estimator 
            hard_simplePrediction = round(mean(value)),
            hard_leaveOneOutPrediction = round(1.0 / mean(weight)),
            leaveOneOutEffectiveSampleSize = (sum(weight))^2 / sum(weight^2)
          )
          
    calibration <- inner_join(dataRaw, presencePosteriors, by = c("locus", "cell"))
    names(calibration)[names(calibration) == 'tipInclusionProbabilities'] <- 'truth'
          
    write.csv(calibration, "«outputFolder.getFileInResultFolder("calibration.csv")»")
    binomial_smooth <- function(...) {
      geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
    }
    p <- ggplot(calibration, aes(leaveOneOutPresence, truth)) +
      binomial_smooth(formula=y~x, alpha=0.2, size=2) +
      geom_abline(intercept = 0) +
      xlab("Leave-one-out predictive probability") + 
      ylab("Binomial smoothed coverage") +
      theme_bw()
    ggsave("«outputFolder.getFileInResultFolder("calibration.pdf")»", p, height = 10, width = 10, limitsize = FALSE)
    
    p <- ggplot(calibration, aes(leaveOneOutPresence)) +
      xlab("Leave-one-out predictive probability") + 
      ylab("Proportion of (cell,locus)") +
      geom_density() +
      theme_bw()
    ggsave("«outputFolder.getFileInResultFolder("uncertainties.pdf")»", p, height = 10, width = 10, limitsize = FALSE)
    
    
    # Hard predictive (i.e. either 0 or 1 - corresponds to 0-1 loss)
        
    hard_predictives <- calibration %>%
      group_by(locus, cell) %>%
      mutate(
        simplePredictive = 1-abs(truth - hard_simplePrediction),
        leaveOneOutPredictive = 1-abs(truth - hard_leaveOneOutPrediction)
      )    
    
    # Predictive (soft call, i.e. in the range [0,1] - corresponds to L2 loss)
        
    obsPredictivesTrace <- raw %>% filter(predicted == 'observation')
    
    predictives <- obsPredictivesTrace %>%
      group_by(locus, cell) %>%
      mutate(weight = 1.0/value) %>%
      summarise(
        simplePredictive = mean(value),
        leaveOneOutPredictive = 1.0 / mean(weight), # harmonic estimator
        leaveOneOutEffectiveSampleSize = (sum(weight))^2 / sum(weight^2)
      )
      
    write.csv(predictives, "«outputFolder.getFileInResultFolder("predictives.csv")»")
    
    # reorg to facilitate plotting
    
    «FOR prediction : #["", "hard_"]»
      
      «prediction»reshaped <- melt(«prediction»predictives, measure.vars = c('leaveOneOutPredictive','simplePredictive'))
      
      «FOR type : #["locus", "cell"]»
      
        # By «type», showing the 2 types of predictives side by side
      
        p <- ggplot(«prediction»reshaped, 
          aes(x = «type», y = value, colour = variable)) +
          geom_boxplot() +
          ylab("Leave-one-out predictive probability") + 
          theme_bw()
        ggsave("«outputFolder.getFileInResultFolder(prediction + "predictives_by_" + type + ".pdf")»", p, height = 10, width = 75, limitsize = FALSE)
        
        «prediction»predictives_by_«type» <- «prediction»predictives %>%
          group_by(«type») %>%
          summarise(
            «prediction»simplePredictiveByLocus = mean(simplePredictive),
            «prediction»leaveOneOutPredictive = mean(leaveOneOutPredictive),
            «prediction»minESS = min(leaveOneOutEffectiveSampleSize)
          )
        write.csv(«prediction»predictives_by_«type», "«outputFolder.getFileInResultFolder(prediction + "predictives_by_" + type + ".csv")»")
        
      «ENDFOR»
      
      # Summary
            
      «prediction»predictives_summary <- «prediction»predictives %>%
        group_by() %>%
        summarise(
          «prediction»simplePredictive = mean(simplePredictive),
          «prediction»leaveOneOutPredictive = mean(leaveOneOutPredictive),
          «prediction»minESS = min(leaveOneOutEffectiveSampleSize)
        )
        
    «ENDFOR»
    
    # Prepare overall summary
        
    baseline <- calibration %>% summarise(value = mean(truth))
    predictives_summary$baseline = max(baseline$value, 1-baseline$value)
    predictives_summary$hard_simplePredictive = hard_predictives_summary$hard_simplePredictive
    predictives_summary$hard_leaveOneOutPredictive = hard_predictives_summary$hard_leaveOneOutPredictive
    predictives_summary$hard_minESS = hard_predictives_summary$hard_minESS
      
    write.csv(predictives_summary, "«outputFolder.getFileInResultFolder("predictives_summary.csv")»")
    
    '''
    callR(outputFolder.getFileInResultFolder(".predictive-script.r"), script)

    // TODO: create matrix with tree
  }
  
  def static void main(String [] args) {
    Experiment::startAutoExit(args)
  }
  
}