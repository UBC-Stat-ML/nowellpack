package corrupt.post

import corrupt.PerfectPhylo
import blang.inits.experiments.tabwriters.TidilySerializable
import blang.inits.experiments.tabwriters.TidySerializer.Context
import corrupt.Locus
import java.util.LinkedHashMap
import corrupt.Cell
import java.util.List
import java.util.ArrayList
import java.util.Collections
import java.util.Random

class PredictiveMonitor implements TidilySerializable {
  val Locus locus
  val BinaryCLMatrix observations
  val PerfectPhylo phylo
  val NoisyBinaryCLMatrix errorModel
  val List<Cell> cells
  
  new (Locus locus, BinaryCLMatrix observations, PerfectPhylo phylo, NoisyBinaryCLMatrix errorModel, double proportion) {
    this.locus = locus
    this.observations = observations
    this.phylo = phylo
    this.errorModel = errorModel
    val _cells = new ArrayList<Cell>(phylo.cells)
    val targetSize = (proportion * phylo.cells.size) as int
    Collections.shuffle(_cells, new Random(1)) // locus.hashCode)) removed locus specific set of cells to ease plotting
    cells = new ArrayList(_cells.subList(0, targetSize))
  }
  
  override serialize(Context context) {
    val tips = phylo.getTips(locus)
    
    val predictiveForObservations = new LinkedHashMap<Cell, Double>
    val predictiveForPresences = new LinkedHashMap<Cell, Double>
    for (cell : cells) {
      val observation = if (observations.get(cell, locus) == 1.0) true else if (observations.get(cell, locus) == 0.0) false else throw new RuntimeException
      val latent = tips.get(cell)
      val predictiveForObservation = p(latent, observation) 
      val predictiveForPresence = p(latent, true)
      
      predictiveForObservations.put(cell, predictiveForObservation)
      predictiveForPresences.put(cell,  predictiveForPresence)
    }
    
    context.recurse(predictiveForObservations, 
      "predicted", "observation") 
    context.recurse(predictiveForPresences, 
      "predicted", "presence") 
  }
  
  def double p(boolean latent, boolean observation) {
    val param = errorModel.parameter(locus)
    if (latent) {
      val fn = errorModel.fn(param)
      if (observation) 1.0 - fn  // 1 -> 1 = 1 - FN
      else             fn        // 1 -> 0 = FN
    } else {
      val fp = errorModel.fp(param)
      if (observation) fp        // 0 -> 1 = FP
      else             1.0 - fp  // 0 -> 0 = 1 - FP
    }
  }
}