package corrupt.post

import blang.runtime.internals.DefaultPostProcessor
import blang.runtime.Runner
import java.io.File
import blang.inits.parsing.ConfigFile
import corrupt.Greedy
import briefj.BriefIO

class CorruptPostProcessor extends DefaultPostProcessor  {
  
  File binaryObservations
  
  var File sampleDir
  override run() {
    // load run
    sampleDir = new File(blangExecutionDirectory.get, Runner::SAMPLES_FOLDER)
    val args = ConfigFile::parse(new File(blangExecutionDirectory.get, DETAILED_ARGUMENT_FILE))
    val binaryObsPath = args.child("model").child("binaryMatrix").argumentValue.get.join
    binaryObservations = new File(binaryObsPath)
    
    // postprocessing steps:
    decode()
    predictiveResults()
    super.run
  }
  
  def decode() {
    val child = results.child("averaging")
    val averager = new AverageCLMatrices => [
      csvFile = new File(sampleDir, "phylo.csv")
      results = child
    ]
    averager.run
    
    val greedyDecoder = new Greedy => [
      tipInclusionProbabilities = ReadOnlyCLMatrix::create(child.getFileInResultFolder(AverageCLMatrices::OUTPUT_NAME).absolutePath) 
    ]
    val consensus = greedyDecoder.infer
    BriefIO::write(results.getFileInResultFolder("consensus.newick"), consensus.toString)
  }
  
  def predictiveResults() {
    val outputFolder = results.child("predictivePlots")
    val predictiveTraces = new File(sampleDir, "predictives.csv") 
    
    val script = '''
    require("ggplot2")
    require("dplyr")
    
    # load data
    
    dataRaw <- read.csv("«binaryObservations»") 
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
            leaveOneOutEffectiveSampleSize = (sum(weight))^2 / sum(weight^2)
          )
          
    calibration <- inner_join(dataRaw, presencePosteriors, by = c("locus", "cell"))
          
    write.csv(calibration, "«outputFolder.getFileInResultFolder("calibration.csv")»")
    binomial_smooth <- function(...) {
      geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
    }
    p <- ggplot(calibration, aes(leaveOneOutPresence, tipInclusionProbabilities)) +
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
    
    
    
    # Predictive stuff
        
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
    
    # By locus, showing the 2 types of predictives side by side
    
    library(reshape2)
    reshaped <- melt(predictives, measure.vars = c('leaveOneOutPredictive','simplePredictive'))
    library(ggplot2)
    
    «FOR type : #["locus", "cell"]»
      p <- ggplot(reshaped, 
        aes(x = «type», y = value, colour = variable)) +
        geom_boxplot() +
        ylab("Leave-one-out predictive probability") + 
        theme_bw()
      ggsave("«outputFolder.getFileInResultFolder("predictives_by_" + type + ".pdf")»", p, height = 10, width = 75, limitsize = FALSE)
      
      predictives_by_«type» <- predictives %>%
        group_by(«type») %>%
        summarise(
          simplePredictiveByLocus = mean(simplePredictive),
          leaveOneOutPredictive = mean(leaveOneOutPredictive),
          minESS = min(leaveOneOutEffectiveSampleSize)
        )
      write.csv(predictives_by_«type», "«outputFolder.getFileInResultFolder("predictives_by_" + type + ".csv")»")
      
      predictives_by_«type» <- predictives %>%
        group_by(«type») %>%
        summarise(
          simplePredictiveByLocus = mean(simplePredictive),
          leaveOneOutPredictive = mean(leaveOneOutPredictive),
          minESS = min(leaveOneOutEffectiveSampleSize)
        )
      write.csv(predictives_by_«type», "«outputFolder.getFileInResultFolder("predictives_by_" + type + ".csv")»")
    «ENDFOR»
    
    # Overall summary
          
    predictives_summary <- predictives %>%
      group_by() %>%
      summarise(
        simplePredictiveByLocus = mean(simplePredictive),
        leaveOneOutPredictive = mean(leaveOneOutPredictive),
        minESS = min(leaveOneOutEffectiveSampleSize)
      )
      
    baseline <- calibration %>% summarise(value = mean(tipInclusionProbabilities))
    predictives_summary$baseline = max(baseline$value, 1-baseline$value)
      
    write.csv(predictives_summary, "«outputFolder.getFileInResultFolder("predictives_summary.csv")»")
    
    '''
    callR(outputFolder.getFileInResultFolder(".predictive-script.r"), script)

    // TODO: create matrix with tree
  }
  
}