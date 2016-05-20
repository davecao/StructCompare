package org.bilab.tools.align;

import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleSubstitutionMatrix;
//import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
//import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.alignment.SubstitutionMatrixHelper;
//import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
//import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.SimpleGapPenalty;

import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;

import java.util.ArrayList;
import java.util.List;

public class SscParameters implements ConfigStrucAligParams {

  protected static final Alignments.PairwiseSequenceAlignerType
       DEFAULT_ALIGN_TYPE = Alignments.PairwiseSequenceAlignerType.LOCAL;
  /*
   *  GLOBAL,              // Needleman-Wunsch/Gotoh
   *  GLOBAL_LINEAR_SPACE, // Guan-Uberbacher
   *  LOCAL,               // Smith-Waterman/Gotoh
   *  LOCAL_LINEAR_SPACE   // Smith-Waterman/Gotoh with smart traceback at each maximum
   */
  protected static final short DEFAULT_GAP_EXTENSION = 1;
  protected static final short DEFAULT_GAP_OPEN = 5;
  protected static final String DEFAULT_ALIGN_ROUTINE="SW_local";
  //protected static final int DEFAULT_ALIGN_ROUTINE = 1;         // "Smith-Waterman/Gotoh"
  protected static final int GLOBAL_ROUTINE = 0;                // "Needleman-Wunsch/Gotoh"
  protected static final int LOCAL_ROUTINE  = 1;                // "Smith-Waterman/Gotoh"
  protected static final int GLOBAL_LINEAR_SPACE_ROUTINE = 2;   // "Guan-Uberbacher"
  protected static final int LOCAL_LINEAR_SPACE_ROUTINE  = 3;   // "Smith-Waterman/Gotoh with smart traceback"
  protected static final String[] ALIGN_ROUTINE = new String[]{"NW_global","SW_local","GU_linear","SW_linear"};
  protected int algoStrategy;// To store the above ROUTINE

  protected SubstitutionMatrix<AminoAcidCompound> substitutionMatrix;

  protected short gapOpen;
  protected short gapExtension;
  protected String algorithmName;
  protected Alignments.PairwiseSequenceAlignerType alignType;
  protected GapPenalty gapPenalty;
  protected String blastAlignFileName = "";
  // Store parsing results from blastp
  //Blast Hsp results
  //protected double hspScore;
  protected double hspEvalue;
  protected double hspBitScore;
  //protected int hspIdentity;
  //protected int hspPositive;
  //protected int hspGaps;
  //protected int hspAlignLen;
  // No Hit found ?
  protected boolean NoHitFound = false;
 
  public SscParameters (){
    reset();
  }
  public SscParameters (String bAlignFileName) {
    blastAlignFileName = bAlignFileName;
    reset();
  }
  // Override the interface ConfigStrucAligParams
  /** Set the parameters to the default.
   *
   */
  public void reset(){
    // SimpleSubstitutionMatrix<AminoAcidCompound>:
    // Blosum 62
    substitutionMatrix = new SimpleSubstitutionMatrix<AminoAcidCompound>();
    gapExtension = DEFAULT_GAP_EXTENSION;
    gapOpen = DEFAULT_GAP_OPEN;
    alignType = DEFAULT_ALIGN_TYPE;
    algorithmName = DEFAULT_ALIGN_ROUTINE;
    algoStrategy = 1;
    // Using the class SimpleGapPenalty.
    gapPenalty = new SimpleGapPenalty();
    gapPenalty.setOpenPenalty(gapOpen);
    gapPenalty.setExtensionPenalty(gapExtension);
  }

  @Override
  public String toString(){
    return "SscParameters [algorithm=" + algorithmName
      + ", gapOpen=" + gapOpen
      + ", gapExtend=" + gapExtension
      + ", BlastAlignmentFile=" + blastAlignFileName
      +"]";
  }

  /** get the list of parameters that the user can change through the user interface.
   *  Parameter names are the same names as the corresponding Get/Set methods.
   *
   * @return list of parameters
   */
  public List<String> getUserConfigParameters(){
    List<String> params = new ArrayList<String>();
    params.add("Algorithm");
    params.add("BlastAlignmentFile");
    //params.add("SubstitutionMatrix");
    params.add("GapOpen");
    params.add("GapExtension");
    //params.add("");
    return params;
  }

  /** The labels to be displayed to the user for each parameter
   *
   * @return list of parameter names
   */
  public List<String> getUserConfigParameterNames(){
    List<String> params = new ArrayList<String>();
    params.add("Algorithm");
    params.add("Using Blast output");
    //params.add("Substitution matrix");
    params.add("Gap Open");
    params.add("Gap Extension");

    return params;
  }

  /** Get the data types of the parameters
   *
   * @return list of parameter classes
   */
  @SuppressWarnings("rawtypes")
  public List<Class> getUserConfigTypes(){
    List<Class> params = new ArrayList<Class>();
    params.add(String[].class);
    params.add(String[].class);
    //params.add(SubstitutionMatrix.class);
    params.add(Short.class);
    params.add(Short.class);

    return params;
  }

  /** The help text for each of these parameters.
   *
   * @return help strings
   */
  public List<String> getUserConfigHelp(){
    List<String> params = new ArrayList<String>();
    params.add("Algorithm:NW_global,SW_local,GU_linear,SW_linear ");
    params.add("Using blast output: the alignment file. The algorithm will be ignored.");
    //params.add("Substitution matrix: scoring matrix used by pairwise sequence alignment");
    params.add("Gap open penalty: used by pairwised sequence alignment");
    params.add("Gap extension penalty: used by pairwised sequence alignment");
    return params;
  }

  /**
   * Return the name of sequence alignment
   *  GLOBAL,              // Needleman-Wunsch/Gotoh
   *  GLOBAL_LINEAR_SPACE, // Guan-Uberbacher
   *  LOCAL,               // Smith-Waterman/Gotoh
   *  LOCAL_LINEAR_SPACE   // Smith-Waterman/Gotoh with smart traceback at each maximum
   * @return algorithm
   */
  public String[] getAlgorithm(){
    return algorithmName.split(":");
  }

  /**
   * Set the name of sequence alignment
   * @param algoName "NW_global","SW_local","GU_linear","SW_linear"
   */
  public void setAlgorithm(String[] aName){

    if ( aName[0].equals("NW_global") ){
      this.algorithmName = aName[0];
      this.alignType = Alignments.PairwiseSequenceAlignerType.GLOBAL;
      this.algoStrategy = 0;
    }else if( aName[0].equals("SW_local") ){
      this.algorithmName = aName[0];
      this.alignType = Alignments.PairwiseSequenceAlignerType.LOCAL;
      this.algoStrategy = 1;
    }else if ( aName[0].equals("GU_linear")  ) {
      this.algorithmName = aName[0];
      this.alignType =  Alignments.PairwiseSequenceAlignerType.GLOBAL_LINEAR_SPACE;
      this.algoStrategy = 2;
    }else if ( aName[0].equals("SW_linear")  ){
      this.algorithmName = aName[0];
      this.alignType = Alignments.PairwiseSequenceAlignerType.LOCAL_LINEAR_SPACE;
      this.algoStrategy = 3;
    }else{
      System.out.println("The input name of the algorithm is unknown: " + aName[0]);
      System.out.println("Use default: "+algorithmName);
    }
  }
  /**
   * Return the blast alignment file name (for GUI use)
   * Default : null
   * @return blastAlignFileName
   */
  public String[] getBlastAlignmentFile(){
    return blastAlignFileName.split("#");
  }

  /**
   * Set the blast alignment file name (for GUI use)
   * @param alignFName
   */
  public void setBlastAlignmentFile(String[] alignFName){
    if (alignFName.length==1){
      blastAlignFileName = alignFName[0];
    }else{
      System.out.println("Unknown file name. This option will be ignored.");
      blastAlignFileName = "";
    }
  }
  /**
   * Return gap extension used in pairwised sequence alignment
   * Default : 0.5
   * @return gapExtension
   */
  public Short getGapExtension(){
    return gapExtension;
  }

  /**
   * Set gap extension used in pairwised sequence alignment
   * Default : 0.5
   * @param gapextension
   */
  public void setGapExtension(Short gapExtension){
    this.gapExtension = gapExtension;
    this.gapPenalty.setExtensionPenalty(gapExtension);
  }

  /**
   * Return gap Open used in pairwised sequence alignment
   * Default : 5.0
   * @return gapOpen
   */
  public Short getGapOpen(){
    return gapOpen;
  }

  /**
   * Set gap extension used in pairwised sequence alignment
   * Default : 5.0
   * @param gapOpen
   */
  public void setGapOpen(Short gapOpen){
    this.gapOpen = gapOpen;
    this.gapPenalty.setOpenPenalty(gapOpen);
  }

  /**
   * Return the substitution matrix
   * Default: BLOSUM62
   * @return substitution
   */
  public SubstitutionMatrix<AminoAcidCompound>
     getSubstitutionMatrix() {
    if ( substitutionMatrix == null){
      this.substitutionMatrix = SubstitutionMatrixHelper.getBlosum62();
    }
    return this.substitutionMatrix;
  }

  /** Sets the substitution matrix.
   * Default: BLOSUM62 matrix
   * @param substitutionMatrix 
   */
  public void setSubstitutionMatrix(
        SubstitutionMatrix<AminoAcidCompound> substitutionMatrix){
    this.substitutionMatrix = substitutionMatrix;
  }

}