package corrupt.post

import bayonet.distributions.Random
import corrupt.PerfectPhylo
import blang.distributions.Normal
import blang.distributions.DiscreteUniform
import blang.distributions.Poisson
import blang.distributions.Binomial

import static blang.types.StaticUtils.*
import bayonet.distributions.Multinomial
import briefj.BriefIO
import java.io.File
import java.util.ArrayList
import corrupt.Locus
import corrupt.Cell
import blang.inits.providers.CoreProviders
import xlinear.Matrix
import blang.distributions.Generators
import java.util.HashSet
import corrupt.distances.LocusEdgeStat
import java.util.Set

import static bayonet.math.NumericalUtils.logAdd


import corrupt.post.SimpleCLMatrix

class CLMatrixUtils {
  
  public static val String CELLS = "cells"
  public static val String LOCI = "loci"
  public static val String TIP_INCL_PRS = "tipInclusionProbabilities"
  public static val String SNV_PROB = 'snv_prob'
  public static val String EVAL_REF = 'evalRef'
  
  static def void toCSV(CellLocusMatrix matrix, File f) {
    toCSV(matrix, f, null)
  } 
  static def void toCSV_snv(CellLocusMatrix matrix, File f, Boolean snv) {
    toCSV_snv(matrix, f, snv, null)
  }
  static def void toCSV(CellLocusMatrix matrix, File f, CellLocusMatrix reference) {
    val out = BriefIO.output(f)
    out.println('''«CELLS»,«LOCI»,«TIP_INCL_PRS»«IF reference !== null»,reference«ENDIF»''')
    for (cell : matrix.cells)
      for (locus : matrix.loci)
        out.println('''«cell»,«locus»,«matrix.get(cell, locus)»«IF reference !== null»,«reference.get(cell, locus)»«ENDIF»''')
    out.close
  }
  
  static def void toCSV_ref(CellLocusMatrix matrix, File f, CellLocusMatrix reference) {
   	  val out = BriefIO.output(f)
	  out.println('''«CELLS»,«LOCI»,«EVAL_REF»«IF reference !== null»,reference«ENDIF»''')
      for (cell : matrix.cells)
        for (locus : matrix.loci){
        	  val locusType = locus.getType()
        	  if (locusType == 'snv')
          	out.println('''«cell»,«locus»,«matrix.get(cell, locus)»«IF reference !== null»,«reference.get(cell, locus)»«ENDIF»''')
      }
      out.close
  } 
  
  static def void toCSV_snv(CellLocusMatrix matrix, File f, Boolean snv, CellLocusMatrix reference) {
    val out = BriefIO.output(f)
    if (!snv){
      out.println('''«CELLS»,«LOCI»,«TIP_INCL_PRS»«IF reference !== null»,reference«ENDIF»''')
      for (cell : matrix.cells)
        for (locus : matrix.loci)
          out.println('''«cell»,«locus»,«matrix.get(cell, locus)»«IF reference !== null»,«reference.get(cell, locus)»«ENDIF»''')
      out.close
    } else {
    	  out.println('''«CELLS»,«LOCI»,«SNV_PROB»«IF reference !== null»,reference«ENDIF»''')
      for (cell : matrix.cells)
        for (locus : matrix.loci){
        	  val locusType = locus.getType()
        	  if (locusType == 'snv')
          	out.println('''«cell»,«locus»,«matrix.get(cell, locus)»«IF reference !== null»,«reference.get(cell, locus)»«ENDIF»''')
      }
      out.close
    }
  }
  
  static def void toCSV_pt(CellLocusMatrix matrix, File f, String pt, CellLocusMatrix reference) {
    val out = BriefIO.output(f)
    out.println('''«CELLS»,«LOCI»,«TIP_INCL_PRS»«IF reference !== null»,reference«ENDIF»''')
    for (cell : matrix.cells)
      for (locus : matrix.loci){
      	if ((locus.getType() == 'cnv') && (locus.getPrintType() == pt)){
            out.println('''«cell»,«locus»,«matrix.get(cell, locus)»«IF reference !== null»,«reference.get(cell, locus)»«ENDIF»''')    
        }else{
      	  if ((locus.getType() == 'snv') && (Float.parseFloat(locus.getPrintType()) <=  Float.parseFloat(pt))){
            out.println('''«cell»,«locus»,«matrix.get(cell, locus)»«IF reference !== null»,«reference.get(cell, locus)»«ENDIF»''')
    		  }         	
        }  
      }
    out.close
    }


  static def int locusEdgeDistance(PerfectPhylo phylo1, PerfectPhylo phylo2) {
    val s1 = LocusEdgeStat::locusEdges(phylo1.collapsedTree)
    val s2 = LocusEdgeStat::locusEdges(phylo2.collapsedTree)
    val union = new HashSet => [
      addAll(s1)
      addAll(s2)
    ]
    val intersection = new HashSet(s1)
    intersection.retainAll(s2)
    union.removeAll(intersection)
    return union.size
  }
  
  static def double distance(PerfectPhylo phylo1, PerfectPhylo phylo2) {
    return distance(fromPhylo(phylo1), fromPhylo(phylo2))
  }
  
  static def double distance(SimpleCLMatrix mtx1, SimpleCLMatrix mtx2) {
    checkCompatible(mtx1, mtx2)
    return distance(mtx1.matrix, mtx2.matrix)
  }
  
  static def double SNVSensitivity(SimpleCLMatrix mtx1, SimpleCLMatrix mtx2, SimpleCLMatrix mtx3) {
    checkCompatible(mtx1, mtx2)
    checkCompatible(mtx1, mtx3)
    return sensitivity(mtx1, round(mtx2), mtx3)
  }
  
   static def double SNVSpecificity(SimpleCLMatrix mtx1, SimpleCLMatrix mtx2, SimpleCLMatrix mtx3) {
    checkCompatible(mtx1, mtx2)
    checkCompatible(mtx1, mtx3)
    return specificity(mtx1, round(mtx2), mtx3)
  }
  
   static def double SNVF1(SimpleCLMatrix mtx1, SimpleCLMatrix mtx2, SimpleCLMatrix mtx3) {
    checkCompatible(mtx1, mtx2)
    checkCompatible(mtx1, mtx3)
    return f1Score(mtx1, round(mtx2), mtx3)
  }
  static def double f1Score(SimpleCLMatrix mtx1, SimpleCLMatrix mtx2, SimpleCLMatrix mtx3) {
	// var double meanF1 = 0.0
	// var int cCount = 0
	var double TP = 0.0
	var double FP = 0.0
	var double TN = 0.0
	var double FN = 0.0	
	// var double pos = 0.0
  	// val double pos =  mtx1.matrix.nonZeroEntries().sum()
  	for (cell : mtx1.cells){
  	  for (locus : mtx1.loci){
  	  	if (mtx3.get(cell,locus) == 1)
  	  	  if (mtx1.get(cell,locus) == 1){
  	        //pos += 1
  	  	    if (mtx2.get(cell,locus) == 1)
  	          TP += 1 
  	        if (mtx2.get(cell,locus) == 0)
  	          FN += 1 
  	  	  }
  	  	  if (mtx1.get(cell,locus) == 0){
  	  	  	if (mtx2.get(cell,locus) == 1)
  	  	  	  FP += 1
  	  	  	if (mtx2.get(cell,locus) == 0)
  	  	  	   TN += 1
  	  	  }
  	  }
  	}  
    return ((2* TP) / (2*TP + FP + FN))
  }
  
   
  static def double sensitivity(SimpleCLMatrix mtx1, SimpleCLMatrix mtx2, SimpleCLMatrix mtx3) {
	//var double meanSensitivity = 0.0
	//var int cCount = 0
	var double TP = 0.0
	var double pos = 0.0
  	// val double pos =  mtx1.matrix.nonZeroEntries().sum()
  	for (cell : mtx1.cells){
  	  //cCount += 1
  	  //TP = 0.0
  	  //pos = 0.0
  	  for (locus : mtx1.loci){
  	  	if (mtx3.get(cell,locus) == 1)
  	  	  if (mtx1.get(cell,locus) == 1){
  	        pos += 1
  	  	    if (mtx2.get(cell,locus) == 1)
  	          TP += 1 
  	  	  }
  	  }
   	  /*if (pos == 0){
  	  	cCount -= 1	
  	  }else{
  	  meanSensitivity += (TP / pos)
  	  //System.out.println('cell sensi = '+  (TP / pos).toString()) 
  	  }	*/
  	}  
  	//System.out.println('final  sensi = '+ (meanSensitivity / cCount).toString()) 
    //return meanSensitivity / cCount
    return TP / pos
  }
  
  static def double specificity(SimpleCLMatrix mtx1, SimpleCLMatrix mtx2, SimpleCLMatrix mtx3) {
  	//var double meanSpecificity = 0.0
	//var int cCount = 0
	var double TN = 0.0
	var double neg = 0.0
  	//val double neg =  mtx1.matrix.nEntries() - mtx1.matrix.nonZeroEntries().sum()
  	for (cell : mtx1.cells){
  	  //cCount += 1	
  	  TN = 0.0
  	  neg = 0.0
  	  for (locus : mtx1.loci){
  	  	if (mtx3.get(cell,locus) == 1)
  	  	  if (mtx1.get(cell,locus) == 0){
  	        neg += 1
  	  	    if (mtx2.get(cell,locus) == 0)
  	          TN += 1 
  	  	  }
  	  }
  	  /*if (neg == 0){
  	  	cCount -= 1	
  	  }else{
  	  meanSpecificity += (TN / neg)	*/
  	  //System.out.println('cell speci = '+  (TN / neg).toString()) 
  	  
  	  //}	
  	}
  	//System.out.println('final speci = '+ (meanSpecificity / cCount).toString()) 
    //return meanSpecificity / cCount
    return TN / neg
  }
  static def double distance(Matrix mtx1, Matrix mtx2) {
    return delta(mtx1, mtx2) / mtx1.nEntries
  }  
  
 
  static def double delta(Matrix mtx1, Matrix mtx2) {
    val diff = mtx1 - mtx2
    return diff.nonZeroEntries().map[Math.abs(it)].sum()
  }
    
  
  static def SimpleCLMatrix fromCSV(File f) {
  	//System.out.println('No boolean-step 1')
  	//System.out.println(f.toString())
    val cells = new ArrayList<Cell>
    val loci = new ArrayList<Locus>
    val values = new ArrayList<Double>
    for (line : BriefIO::readLines(f).indexCSV) {
      if (!line.containsKey(CELLS) || !line.containsKey(LOCI) || !line.containsKey(TIP_INCL_PRS))
        throw new RuntimeException("The following column headers should be present in " + f.name + ": " + 
          CELLS + ", " + LOCI + ", " + TIP_INCL_PRS + "\nFound: " + line.keySet.join(", ")
        )
        
      
      cells.add(new Cell(line.get(CELLS)))
      loci.add(new Locus(line.get(LOCI)))
      values.add(CoreProviders::parse_double(line.get(TIP_INCL_PRS)))
    }

    
    val result = new SimpleCLMatrix(cells, loci)
    for (var int i = 0; i < cells.size; i++) 
      result.set(cells.get(i), loci.get(i), values.get(i))
    return result
  }
  
  
  static def SimpleCLMatrix fromCSV(File f, Boolean snv) {
  	//System.out.println('from boolean')
    val cells = new ArrayList<Cell>
    val loci = new ArrayList<Locus>
    val values = new ArrayList<Double>
    for (line : BriefIO::readLines(f).indexCSV) {
    	  if (!snv){
        if (!line.containsKey(CELLS) || !line.containsKey(LOCI) || !line.containsKey(TIP_INCL_PRS))
          throw new RuntimeException("The following column headers should be present in " + f.name + ": " + 
            CELLS + ", " + LOCI + ", " + TIP_INCL_PRS + "\nFound: " + line.keySet.join(", ")
          )
        cells.add(new Cell(line.get(CELLS)))
        loci.add(new Locus(line.get(LOCI)))
        values.add(CoreProviders::parse_double(line.get(TIP_INCL_PRS)))
      } else {
        if (!line.containsKey(CELLS) || !line.containsKey(LOCI))
          throw new RuntimeException("The following column headers should be present in " + f.name + ": " + 
            CELLS + ", " + LOCI + "\nFound: " + line.keySet.join(", ")
          )
        if (!line.containsKey(SNV_PROB) && !line.containsKey(EVAL_REF))
          throw new RuntimeException("The following column headers should be present in " + f.name + ": " + 
            EVAL_REF + " or " + SNV_PROB + "\nFound: " + line.keySet.join(", ")
          )
        cells.add(new Cell(line.get(CELLS)))
        loci.add(new Locus(line.get(LOCI)))
        
        if (line.containsKey(SNV_PROB))
        	  values.add(CoreProviders::parse_double(line.get(SNV_PROB)))
        	if (line.containsKey(EVAL_REF))
        	  values.add(CoreProviders::parse_double(line.get(EVAL_REF)))
      }
    }
    val result = new SimpleCLMatrix(cells, loci)
    for (var int i = 0; i < cells.size; i++) 
      result.set(cells.get(i), loci.get(i), values.get(i))
    return result
  }
  
  static def SimpleCLMatrix fromPhylo(PerfectPhylo phylo) {
    val result = new SimpleCLMatrix(phylo.cells, phylo.loci)
    result += phylo 
    return result
  }
  static def SyntheticReferenceMatrix(Random rand, PerfectPhylo phylo){
  	SyntheticReferenceMatrix(rand, phylo, phylo.loci)
  }
  static def SyntheticReferenceMatrix(Random rand, PerfectPhylo phylo, Set<Locus> loci){
  	val referencePrs = new SimpleCLMatrix(phylo.cells, loci)
  	for (locus : loci){
  	  val tips = phylo.getTips(locus)
  	  for (cell:phylo.cells){
  	    var mutStatus = 0
  	    var ref = tips.get(cell)
  	    if (ref)
  	      mutStatus =1	
  	    referencePrs.set(cell, locus, mutStatus)
  	  }
  	}
  	return ReadOnlyCLMatrix.readOnly(referencePrs)
  }
  /*static def SimpleCLMatrix fromPhylo(PerfectPhylo phylo, String pt) {
    val result = new SimpleCLMatrix(phylo.cells, phylo.loci, pt)
    result.fillPhyloMatrix(phylo, pt) 
    return result
  }*/
  
  static def ReadOnlyCLMatrix syntheticInclusionPrs(Random rand, PerfectPhylo phylo, double stdDev, double seqError, double seqDO, double coverage, int cMax, File f) {
    syntheticInclusionPrs(rand, phylo, stdDev, phylo.loci, seqError, seqDO, coverage, cMax, f)
  }
  static def ReadOnlyCLMatrix syntheticInclusionPrs(Random rand, PerfectPhylo phylo, double stdDev, Set<Locus> loci, double seqError, double seqDO, double coverage, int cMax, File f) {
    val syntheticInclusionPrs = new SimpleCLMatrix(phylo.cells, loci)
    val evalRef = new SimpleCLMatrix(phylo.cells, loci)
    
    val inclNormal = Normal::distribution(fixedReal(0.0), fixedReal(stdDev * stdDev))
    val exclNormal = Normal::distribution(fixedReal(1.0), fixedReal(stdDev * stdDev))
    for (locus : loci) {
    	  val type = locus.getType()
    	  val tips = phylo.getTips(locus)
    	  if (type == 'cnv'){
        for (cell : phylo.cells) {
          val dist = if (tips.get(cell)) inclNormal else exclNormal
          val observation = dist.sample(rand)
          val prs = newDoubleArrayOfSize(2)
          prs.set(0, exclNormal.logDensity(observation))
          prs.set(1, inclNormal.logDensity(observation))
          Multinomial::expNormalize(prs)
          syntheticInclusionPrs.set(cell, locus, prs.get(1)) 
          evalRef.set(cell, locus, 0) 
        } 
      }
      if (type == 'snv'){
        val double one = 1
      	for (cell:phylo.cells){
      	  evalRef.set(cell, locus, 1)
      	  var cn = DiscreteUniform::distribution(fixedInt(0), fixedInt(cMax))
      	  var c = cn.sample(rand)
      	  var prs = newDoubleArrayOfSize(2)
      	  var dd = coverage 
      	  if (c>0){
      	    dd = dd*c/2
      	  }
          var depth = Poisson::distribution(fixedReal(dd))
      	  var d = depth.sample(rand)      		
          var b = 0
          if (d == 0){ 
          	evalRef.set(cell, locus, 0) 		
      		prs.set(1, Math.log(0.5))
			prs.set(0, Math.log(0.5))
		  } else {	
      		if (c == 0){
      		  evalRef.set(cell, locus, 0)	
      		  prs.set(1, Math.log(seqError))
      		  prs.set(0, Math.log(1- seqError))
      		} else { 
      		  if (!tips.get(cell)){  
      		    /*var binomSeqError = Binomial::distribution(fixedInt(d), fixedReal(seqError))
      			b = binomSeqError.sample(rand)    		  
      		    prs.set(1, binomSeqError.logDensity(b))
      		    prs.set(0, Math.log(1 - Math.exp(binomSeqError.logDensity(b))))*/
      		    prs.set(1, Math.log(seqError))
      		    prs.set(0, Math.log(1- seqError))
      		  }else{     		
      		    /*var DODist = Bernoulli::distribution(fixedReal(seqDO))*/
      		    var DODist = Binomial::distribution(fixedInt(1), fixedReal(seqDO))
      		    var DO = DODist.sample(rand)
      		    /*System.out.println ("DO  equals " + DO.toString)*/
      		    if (DO.toString == '1'){
      		      /*System.out.println ("DO = 1.000000")*/
      		      b = 0
      		    } else {
      		    	  /*System.out.println ("Not DO == 0.0000" )*/
      		      var g = 1000000
      			  if (c.toString == '1'){
      			    g = 1
      			  } else{     		      
      		        var gDist = DiscreteUniform::distribution(fixedInt(1), fixedInt(c))	
      		        g = gDist.sample(rand)
      		      }
      		      /*System.out.println ("g = " + g.toString)*/
      		      var xii = (one * g)   / c
      		      var binomB = Binomial::distribution(fixedInt(d), fixedReal(xii)) 
      		      b = binomB.sample(rand)
      		      /*System.out.println ("b = " + b.toString )*/
                }
        		  	var double pg = (1 - seqDO) / (c)
				var binom = Binomial::distribution(fixedInt(d), fixedReal(seqError))
           		var incpr = Math.log(seqDO) + binom.logDensity(b)
      		    binom = Binomial::distribution(fixedInt(d), fixedReal(1-seqError)) 
      		    incpr = logAdd(incpr, Math.log(pg) + binom.logDensity(b))
      		    /*System.out.println ("pr1 = " + Math.exp(binom.logDensity(b) + Math.log(pg)).toString )*/
      		    for (var int i=1; i<c; i++){ 
      		    	  var xi = (one * i)/c
      		      binom = Binomial::distribution(fixedInt(d), fixedReal(xi)) 
      		  	  incpr = logAdd(incpr, Math.log(pg) + binom.logDensity(b))	
      		    }	
      		    prs.set(1, incpr)
      		    prs.set(0, Math.log(1 - Math.exp(incpr)))
      		  }
      		}
          }
          syntheticInclusionPrs.set(cell, locus,  Math.exp(prs.get(1)))
        }
      }
    }  
    ReadOnlyCLMatrix.readOnly(evalRef).toCSV_ref(f, CLMatrixUtils::fromPhylo(phylo)) 
    return ReadOnlyCLMatrix.readOnly(syntheticInclusionPrs)
  }
  
  static def BinaryCLMatrix syntheticPerturbedBinaryMatrix(Random rand, PerfectPhylo phylo, double fpRate, double fnRate) {
    syntheticPerturbedBinaryMatrix(rand, phylo, fpRate, fnRate, phylo.loci)
  }
  static def BinaryCLMatrix syntheticPerturbedBinaryMatrix(Random rand, PerfectPhylo phylo, double fpRate, double fnRate, Set<Locus> loci) {
    val syntheticInclusionPrs = new SimpleCLMatrix(phylo.cells, loci)
    for (locus : loci) {
      val tips = phylo.getTips(locus)
      for (cell : phylo.cells) {
        var indic = tips.get(cell)
        if (indic) {
          if (Generators::bernoulli(rand, fnRate))
            indic = false
        } else {
          if (Generators::bernoulli(rand, fpRate))
            indic = true
        }
        syntheticInclusionPrs.set(cell, locus, if (indic) 1.0 else 0.0) 
      }  
    }
    return BinaryCLMatrix::binary(syntheticInclusionPrs)
  } 
    
  static def ReadOnlyCLMatrix syntheticPerturbedBinaryMatrix(Random rand, PerfectPhylo phylo, double fpRate, double fnRate, double seqError, double seqDO, double coverage, int cMax, File f) {
    syntheticPerturbedBinaryMatrix(rand, phylo, fpRate, fnRate, phylo.loci, seqError, seqDO, coverage, cMax, f)
  }
  static def ReadOnlyCLMatrix syntheticPerturbedBinaryMatrix(Random rand, PerfectPhylo phylo, double fpRate, double fnRate, Set<Locus> loci, double seqError, double seqDO, double coverage, int cMax, File f) {
    val syntheticInclusionPrs = new SimpleCLMatrix(phylo.cells, loci)    
    val evalRef = new SimpleCLMatrix(phylo.cells, loci)
    for (locus : loci) {
      val tips = phylo.getTips(locus)
      val type = locus.getType()
      if (type == 'cnv'){
        for (cell : phylo.cells) {
          var indic = tips.get(cell)
          if (indic) {
            if (Generators::bernoulli(rand, fnRate))
              indic = false
          } else {
            if (Generators::bernoulli(rand, fpRate))
              indic = true
          }
          syntheticInclusionPrs.set(cell, locus, if (indic) 1.0 else 0.0) 
          evalRef.set(cell, locus, 0) 
          
        }
      }
      if (type == 'snv'){
        val double one = 1
      	for (cell:phylo.cells){
      	  //System.out.println("new cell")
      	  evalRef.set(cell, locus, 1) 	
      	  var cn = DiscreteUniform::distribution(fixedInt(1), fixedInt(cMax))
      	  //System.out.println("cn = " + cn.toString())
      	  var c = cn.sample(rand)
      	  var prs = newDoubleArrayOfSize(2)
      	  var dd = coverage
      	  if (c>0){
      	    dd = (coverage/2.0)*c
      	  }
      	  //System.out.println("coverage = " + dd.toString())
      	  //System.out.println("c = " + c.toString())
      	  
          var depth = Poisson::distribution(fixedReal(dd))
      	  var d = depth.sample(rand)
      	  //System.out.println("d = " + d.toString())
          var b = 0
          if (d == 0){
          	// evalRef.set(cell, locus, 0) 
        		prs.set(1, Math.log(0.5))
			prs.set(0, Math.log(0.5))
		  } else {	
      		if (c == 0){
              // ievalRef.set(cell, locus, 0) 
      		  prs.set(1, Math.log(seqError))
      		  prs.set(0, Math.log(1- seqError))
      		} else { 
      		  if (!tips.get(cell)){ 	  
      		    prs.set(1, Math.log(seqError))
      		    prs.set(0, Math.log(1- seqError))
      		  } else { 
      		    var DODist = Binomial::distribution(fixedInt(1), fixedReal(seqDO))
      		    var DO = DODist.sample(rand)
      		    if (DO.toString == '1'){
      		      b = 0
      		    } else {
      		      var g = 1000000
      			  if (c.toString == '1'){
      			    g = 1
      			  } else{     		      
      		        var gDist = DiscreteUniform::distribution(fixedInt(1), fixedInt(c))	
      		        g = gDist.sample(rand)
      		      }
      		      var xii = (one * g)   / c
      		      var binomB = Binomial::distribution(fixedInt(d), fixedReal(xii)) 
      		      b = binomB.sample(rand)
      		      //System.out.println("b = " + b.toString())
                }
				var binom = Binomial::distribution(fixedInt(d), fixedReal(seqError))
				var incpr0 = binom.logDensity(b)
           		var incpr_DO0 = Math.log(seqDO) + binom.logDensity(b)
           		//System.out.println("tip without DO = " + Math.exp(incpr_DO0).toString())
           		
          	    var double pg = (1 - seqDO) / (c)
      		    binom = Binomial::distribution(fixedInt(d), fixedReal(1-seqError)) 
      		    var incpr1 = logAdd(incpr_DO0, Math.log(pg) + binom.logDensity(b) )		
      		    // var temp_pr =  binom.logDensity(b)
      		    //System.out.println("pr for xi=1 " + Math.exp(temp_pr).toString())
      		    for (var int i=1; i<c; i++){ 
      		    	  var xi = (one * i)/c
      		    	  //System.out.println("xi = " + xi.toString())	  
      		      binom = Binomial::distribution(fixedInt(d), fixedReal(xi)) 
      		      incpr1 = logAdd(incpr1, Math.log(pg) + binom.logDensity(b))	
      		      //System.out.println("prob for (i = " + i.toString() + ") = " +  Math.exp( binom.logDensity(b)).toString())
      		      //if (binom.logDensity(b)>temp_pr){
      		      	//temp_pr = binom.logDensity(b)
      		      	//System.out.println("choosed "+ i.toString())
      		      //}      		  	  
      		    }	
      		    //var incpr_DO1 = logAdd(Math.log(1- seqDO) , temp_pr)	
      		    //System.out.println("tip with DO = " + Math.exp(incpr_DO1).toString())
      		    
      		    
      		    //var incpr =  logAdd(incpr_DO0, incpr_DO1)
      		    prs.set(1, incpr1 - logAdd(incpr1 , incpr0) )
      		    prs.set(0, incpr0 - logAdd(incpr1 , incpr0))
      		    
      		    
      		    //System.out.println("tip = " + Math.exp(incpr).toString())
      		    //prs.set(0, Math.log(1 - Math.exp(incpr)))
      		   /* if (b>0){
      		    		prs.set(0, Math.log(seqError))
      		    		prs.set(1, Math.log(1- seqError))
      		    }else{
      		    		prs.set(1, Math.log(seqError))
      		    		prs.set(0, Math.log(1- seqError))    
      		    } */
      		    //System.out.println("final incrp1" +   Math.exp(prs.get(1)).toString())		    	
      		    
      		  }
      		}
          }
          syntheticInclusionPrs.set(cell, locus, Math.exp(prs.get(1)))
        }
      }
    }
    ReadOnlyCLMatrix.readOnly(evalRef).toCSV_ref(f, CLMatrixUtils::fromPhylo(phylo)) 
    return ReadOnlyCLMatrix.readOnly(syntheticInclusionPrs)
  } 
  
  static def void checkCompatible(SimpleCLMatrix cl1, SimpleCLMatrix cl2) {
    if (cl1.cellsIdx != cl2.cellsIdx || cl1.lociIdx != cl2.lociIdx){
    	  throw new RuntimeException
    }
  }
  
  static def SimpleCLMatrix round(SimpleCLMatrix clm ) {
  	var binaryMat = new SimpleCLMatrix(clm.cells, clm.loci)
    for (cell : clm.cells) 
      for (locus : clm.loci) 
        binaryMat.set(cell, locus, Math.round(clm.get(cell, locus))) 
    return binaryMat
  }
}