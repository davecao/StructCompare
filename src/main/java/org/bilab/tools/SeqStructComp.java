package org.bilab.tools;

//import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
//import org.biojava3.alignment.template.SequencePair;
//import org.biojava3.alignment.template.SubstitutionMatrix;
//import org.biojava3.core.sequence.ProteinSequence;
//import org.biojava3.core.sequence.compound.AminoAcidCompound;
//import org.biojava3.core.sequence.io.FastaReaderHelper;

import org.biojava.bio.structure.*;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.AbstractStructureAlignment;
import org.biojava.bio.structure.align.StructureAlignment;
//import org.biojava.bio.structure.io.PDBFileReader;
//import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;
//import org.biojava.bio.structure.jama.Matrix;

import org.bilab.tools.align.SscParameters;
import org.bilab.tools.align.SSCCalculator;

public final class SeqStructComp
  extends AbstractStructureAlignment
       implements StructureAlignment {

  public static String algorithmName = "jSeqBase";
  public static final String version = "0.1";
  public static String newline = System.getProperty("line.separator");

  protected SscParameters params;
  protected SSCCalculator calculator;
  private Atom[] ca2clone;
  /**
   * Construct a SeqStructComp object
   */
  public SeqStructComp(){
    params = new SscParameters();
  }
    /*  static factory method vs constructor
  public void setArg1(String mol){
    mol1 = mol;
  }
  public void setArg2(String chain){
    chain1 = chain;
  }

  public static SeqStructComp create(String mola,String chaina,String molb,String chainb){
    SeqStructComp c = new SeqStructComp();
    c.setArg1(mol1);
    c.setArg2(chain1);
    return c;
  }
  */

  /**
   * Return the algorithm name
   * @return algorithmName
   */
  public String getAlgorithmName(){
    return algorithmName;
  }

  /**
   * Return the version info
   * @return version
   */
  public String getVersion(){
    return version;
  }

  /** Returns some documentation on the command line arguments for this algorithm.
   *
   * @return the help string
   */
  public String printHelp(){
    StringBuffer buf = new StringBuffer();
    buf.append("-------------------").append(newline);
    buf.append("jSeqBase v." + version  + " help: " + newline);
    buf.append("-------------------").append(newline);
    return buf.toString();
  }

  /** Run an alignment while specifying the atoms to be aligned. Will used default parameters for the algorithm.
   *
   * @param ca1
   * @param ca2
   * @return the afpChain object that contains the alignment.
   * @throws StructureException
   */

  public AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException{
     if (params == null)
            params = new SscParameters();
     return align(ca1,ca2,params);
  }

  /** run an alignment and also send a bean containing the parameters.
   *
   * @param ca1
   * @param ca2
   * @param params
   * @return the afpChain object that contains the alignment.
   * @throws StructureException
   */

  public AFPChain align(Atom[] ca1, Atom[] ca2, Object param) throws StructureException{
    if ( ! (param instanceof SscParameters))
          throw new IllegalArgumentException("We needs an object of call SscParameters as argument.");
    params = (SscParameters) param;
    ca2clone = new Atom[ca2.length];
    int pos = 0;
    for (Atom a : ca2){
      Group g = (Group)a.getGroup().clone();// works because each group has only a CA
      ca2clone[pos] = g.getAtom(StructureTools.caAtomName);
      pos++;
    }
    calculator = new SSCCalculator(params);
    // Step 4. build a structural alignment.
    AFPChain afpChain = new AFPChain();
    afpChain = calculator.align(ca1, ca2clone);
    //calculator.traceFragmentMatrix( afpChain,ca1, ca2clone);
    //calculator.nextStep( afpChain,ca1, ca2clone);

    afpChain.setAlgorithmName(getAlgorithmName());
    afpChain.setVersion(version);
    // Step 5. Set the distance matrix.
    //double[][] m = calculator.iniSumOfDistances(ca1.length,ca2.length,ca1,ca2clone); // Not completed.
    //afpChain.setDistanceMatrix(new Matrix(m));
    //afpChain.setSequentialAlignment(true);
    return afpChain;
  }

  /**
   * Set the default parameters for this algorithm to use
   *
   * @param parameters
   */

  public void setParameters(ConfigStrucAligParams parameters){
    params = (SscParameters)parameters;
  }

 /** Return the paramers for this algorithm.
  *
  * @return The returned object will be a Java bean.
  */

  public ConfigStrucAligParams getParameters(){
    return params;
  }

}