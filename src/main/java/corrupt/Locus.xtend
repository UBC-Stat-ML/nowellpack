package corrupt

import blang.inits.DesignatedConstructor
import blang.inits.Input

class Locus extends TreeNode {
  String type
  String printType
  public static val PREFIX = "locus_"
  
  @DesignatedConstructor new(@Input String id, String t, String pt) { 
  	super(id.replaceFirst(PREFIX,  "")) 
    type = t
    printType = pt
  }
  
  @DesignatedConstructor new(@Input String id) { 
  	super(id.replaceFirst(PREFIX,  ""))
    type = 'cnv'
    printType = '0'
  }
  override toString() { PREFIX + id }
  
  
  def String getPrintType(){
  	return printType
  }
  def void setPrintType(String pt){
  	this.printType = pt
  }
  
  def String getType(){
  	return type
  }
  def void setType(String t){
  	this.type = t
  }
}

/* 
class CNVLocus extends Locus {
  public static val PREFIX = "locus_"
  @DesignatedConstructor new(@Input String id) { super(id.replaceFirst(PREFIX,  "")) }
  override toString() { PREFIX + id }
}

class SNVLocus extends Locus {
  public static val PREFIX = "locus_"
  @DesignatedConstructor new(@Input String id) { super(id.replaceFirst(PREFIX,  "")) }
  override toString() { PREFIX + id }
}
*/
