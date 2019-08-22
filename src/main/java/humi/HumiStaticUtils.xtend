package humi

import blang.inits.Creator
import blang.inits.Creators
import blang.inits.providers.CoreProviders
import blang.io.Parsers
import blang.io.internals.GlobalDataSourceStore
import blang.inits.experiments.ParsingConfigs
import blang.inits.experiments.ExperimentResults
import briefj.BriefIO
import humi.freq.DeltaMethod.Columns
import binc.Command

class HumiStaticUtils {
  
  def static plotIntervals(ExperimentResults results, HumiData data, String rCmd, boolean expSpecific, String caption) {
    val plotResults = results.child("plots")
    val scriptFile = plotResults.getFileInResultFolder("script.r")
    BriefIO::write(scriptFile, '''
      require("ggplot2")
      data <- read.csv("«results.getFileInResultFolder("estimates.csv").absolutePath»")
      cols = rainbow(200, s=.6, v=.9)[sample(1:200,200)]
      p <- ggplot(data, aes(x = factor(«data.genes.name»), y = «Columns::logRatio», colour = factor(«data.targets.name»))) + 
        coord_flip() + 
        geom_errorbar(aes(ymin=«Columns::logRatioLeftBound», ymax=«Columns::logRatioRightBound»)) +
        geom_point() + 
        «IF expSpecific»facet_grid(. ~ «data.experiments.name») + «ENDIF»
        theme_bw() + 
        xlab("Gene") + 
        ylab("log(ratio)") + 
        ggtitle("Ratio of clone sizes relative to controls", subtitle = "«caption»") + 
        scale_colour_manual(values=cols) + 
        geom_hline(yintercept=0) + 
        theme(legend.position="none") 
      ggsave("«plotResults.getFileInResultFolder("intevals.pdf").absolutePath»", height = 10)
    ''')
    Command.call(Command.cmd(rCmd).appendArg(scriptFile.getAbsolutePath()))
  }

  // TODO: move to Blang SDK
  def static parsingConfigs() {
    val Creator creator = Creators::empty()
    creator.addFactories(CoreProviders)
    creator.addFactories(Parsers)
    val GlobalDataSourceStore globalDS = new GlobalDataSourceStore
    creator.addGlobal(GlobalDataSourceStore, globalDS)
    val ParsingConfigs parsingConfigs = new ParsingConfigs
    parsingConfigs.setCreator(creator) 
    return parsingConfigs
  }
}