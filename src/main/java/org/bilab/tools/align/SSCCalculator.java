package org.bilab.tools.align;

import java.lang.StringBuilder;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.SVDSuperimposer;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
//import org.biojava.bio.structure.align.model.AFP;
import org.biojava.bio.structure.align.model.AFPChain;
//import org.biojava.bio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.bio.structure.jama.Matrix;
//import org.biojava.bio.structure.align.ce.MatrixListener;

import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.ce.CeParameters;
//import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;
//import org.biojava.bio.structure.align.ce.UserArgumentProcessor;

import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.PairwiseSequenceAligner;
//import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
//import org.biojava3.alignment.aaindex.ScaledSubstitutionMatrix;
//import org.biojava3.alignment.template.SubstitutionMatrix;

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.template.Compound;

import org.bilab.tools.align.SscParameters;
import org.bilab.tools.align.BlastXMLParser;

public class SSCCalculator {

  double[][] dist1;
  double[][] dist2;
  long timeStart;
  long timeEnd;

  //  private int nAtom;

  //  private int[] align_se1;
  //  private int[] align_se2;
  //  List<MatrixListener> matrixListeners;

  //  protected double[][]	mat;
  protected Matrix rotMatrix;
  protected Atom  tranMatrix;
  protected SscParameters params;

  public SSCCalculator(SscParameters params){
    timeStart = System.currentTimeMillis();
    //    dist1= new double[0][0];
    //    dist2= new double[0][0];
    this.params = params;
    //    matrixListeners = new ArrayList<MatrixListener>();
  }

  public AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException {
    if ( params == null){
      params = new SscParameters();
    }
    if (params.blastAlignFileName.isEmpty()){
    	return align(ca1,ca2,params);
    }
    return align(ca1,ca2,params.blastAlignFileName);
  }
  /**
   * Using specified blast output to do structural alignment
   *
   * @param ca1
   * @param ca2
   * @param blastAlignFName the file of blast output (xml)
 * @throws StructureException 
   */
  public AFPChain align(Atom[] ca1, Atom[] ca2, String blastAlignFName) throws StructureException {
    AFPChain afpChain = new AFPChain();
    // Step 1. get the pair-wised alignment
    /* Call new object
    BlastXMLParser<ProteinSequence, AminoAcidCompound> bparse = new BlastXMLParser<ProteinSequence, AminoAcidCompound>(); 
    bparse.xmlFile = blastAlignFName;
    PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> seqAlignment 
        = bparse.parse();
    */
    // Step 1. extract  sequences from structures.
    String seq1 = StructureTools.convertAtomsToSeq(ca1);
    String seq2 = StructureTools.convertAtomsToSeq(ca2);

    ProteinSequence orignal1 = new ProteinSequence(seq1);
    ProteinSequence orignal2 = new ProteinSequence(seq2);
    
    // Call from static factory method
    PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> seqAlignment
           = BlastXMLParser.parseBlastXML(blastAlignFName,orignal1,orignal2);
    // Step 2. find corresponding pairs.
    SequencePair<ProteinSequence, AminoAcidCompound> pair = seqAlignment.getPair();
   
    List<Atom> nCa1 = new ArrayList<Atom>();
    List<Atom> nCa2 = new ArrayList<Atom>();
    int nPairs = pair.getLength();
    int nIdentical = 0;
    int nSimilar   = 0;
    int nGaps      = 0;
    int numNonGapCols = 0;
    StringBuilder alnQuery  = new StringBuilder();
    StringBuilder alnTarget = new StringBuilder();
    StringBuilder alnSymbol = new StringBuilder();

    Compound gapSymbol =  AminoAcidCompoundSet.getAminoAcidCompoundSet().getCompoundForString("-");
    int[] align_se1 = new int[nPairs+1];
    int[] align_se2 = new int[nPairs+1];

    for (int i=1; i<=nPairs; i++){
      int currPos = i-1;
      Compound s1 = pair.getCompoundAt(1,i);
      Compound s2 = pair.getCompoundAt(2,i);
      //Position in original ca1 or ca2
      int Pos1 = pair.getIndexInQueryAt(i)-1;
      int Pos2 = pair.getIndexInTargetAt(i)-1;

      if ( !s1.equals(gapSymbol) && !s2.equals(gapSymbol) ){
        // normal amino acid pairs in the column
        numNonGapCols++;

        nCa1.add(ca1[Pos1]);
        nCa2.add(ca2[Pos2]);

        alnQuery.append(s1.getShortName());
        alnTarget.append(s2.getShortName());

        if(s1==s2){//An identical pair
          nIdentical++;
          nSimilar++;
          alnSymbol.append("|");
        }else{     //A similar pair
          nSimilar++;
          alnSymbol.append(":");
        }
        //Record the position in the array of ca1 or ca2
        align_se1[currPos] = Pos1;
        align_se2[currPos] = Pos2;

      }else if( s1.equals(gapSymbol) && s2.equals(gapSymbol) ){
        // There are all gaps in the column
        nGaps = nGaps+2;
        alnQuery.append("-");
        alnTarget.append(" ");
        alnSymbol.append(" ");
        align_se1[currPos] = -1;
        align_se2[currPos] = -1;
      }else {
        // There is a gap at least in the column
        nGaps++;
        alnSymbol.append(" ");

        align_se1[currPos] = -1;
        align_se2[currPos] = -1;

        if ( s1.equals(gapSymbol) ){
          alnQuery.append("-");
          alnTarget.append(s2.getShortName());

          align_se2[currPos] = Pos2;

        }else{
          alnQuery.append(s1.getShortName());
          alnTarget.append("-");

          align_se1[currPos] = Pos1;
        }
      }
    }
    // Print out the alignment
    //System.out.println(alnQuery.toString());
    //System.out.println(alnSymbol.toString());
    //System.out.println(alnTarget.toString());
    printPSAinfo(seqAlignment);

    afpChain.setAlnLength(nPairs);
    afpChain.setGapLen(nGaps);
    afpChain.setAlnseq1(alnQuery.toString().toCharArray());
    afpChain.setAlnseq2(alnTarget.toString().toCharArray());
    afpChain.setAlnsymb(alnSymbol.toString().toCharArray());

    afpChain.setIdentity(nIdentical*1.0/numNonGapCols);
    afpChain.setSimilarity(nSimilar*1.0/numNonGapCols);
    afpChain.setOptLen(new int[]{numNonGapCols});
    afpChain.setAlignScore(seqAlignment.getScore());
    afpChain.setOptLength(numNonGapCols);
//    afpChain.setProbability();
    // original length of ca1 and ca2
    afpChain.setCa1Length(ca1.length);
    afpChain.setCa2Length(ca2.length);
    //afpChain.setCa1Length(nCa1.size());
    //afpChain.setCa2Length(nCa2.size());

    //Less than 4 atoms
    if (alnQuery.length() < 4){
      return afpChain;
    }

    // do a 3D alignment...
    Atom[] qCa1 = nCa1.toArray(new Atom[0]);
    Atom[] tCa2 = nCa2.toArray(new Atom[0]);

    //System.out.println("Query  atoms: "+qCa1.length);
    //System.out.println("Target atoms: "+tCa2.length);

    // SVD to get Rotation coords and Shift coords
    SVDSuperimposer svds = new SVDSuperimposer(qCa1, tCa2);
    this.rotMatrix = svds.getRotation();
    this.tranMatrix = svds.getTranslation();

    // now we have all the info to perform the rotations ...
    for (Atom X : tCa2 ){
      Calc.rotate(X,this.rotMatrix);
      // shift Ca2 onto Ca1 ...
      Calc.shift(X,this.tranMatrix);
    }
    // Assert
    assert(qCa1.length == tCa2.length);
    double rmsd = SVDSuperimposer.getRMS(qCa1, tCa2);
    // Configure the  afpChain
    afpChain.setBlockRmsd(new double[]{rmsd});
    afpChain.setOptRmsd(new double[]{rmsd});
    afpChain.setTotalRmsdOpt(rmsd);
    afpChain.setChainRmsd(rmsd);

    //Using ce: refer to the sample SmithWaterman3Daligner.java
    CeParameters params = new CeParameters();
    CECalculator cecalc = new CECalculator(params);
    cecalc.setnAtom(numNonGapCols);
    //align_se1 and align_se2 are necessary for convertAfpChain(...)
    cecalc.setAlign_se1(align_se1);
    cecalc.setAlign_se2(align_se2);
    //cecalc.setMatMatrix(this.rotMatrix.getArray());
    cecalc.setLcmp(nPairs);
    // Create corresponding pairs from Structure info
    // utilize the routine in ce
    cecalc.convertAfpChain(afpChain, ca1, ca2);
    // Elasped time
    this.timeEnd = System.currentTimeMillis();
    afpChain.setIoTime(this.timeEnd-this.timeStart);
    // save the rotation matrix and the shift vector
    afpChain.setBlockRotationMatrix(new Matrix[]{this.rotMatrix});
    afpChain.setBlockShiftVector(new Atom[]{this.tranMatrix});

    return afpChain;
  }

  /**
   * Using pairwise alignment algorithm in Biojava
   * @param ca1
   * @param ca2
   * @param parameters SscParameters
   */
  public AFPChain align(Atom[] ca1, Atom[] ca2, Object parameters) throws StructureException{
    if ( parameters == null) {
           throw new IllegalArgumentException("Got null instead of SscParameters!");
    }
    if ( ! (parameters instanceof SscParameters)){
              throw new IllegalArgumentException("provided parameter object is not of type SscParameters,but "
                      +parameters.getClass().getName());
    }
    params = (SscParameters) parameters;
    AFPChain afpChain = new AFPChain();
    try {
      // Step 1. extract  sequences from structures.
      String seq1 = StructureTools.convertAtomsToSeq(ca1);
      String seq2 = StructureTools.convertAtomsToSeq(ca2);

      ProteinSequence protein1 = new ProteinSequence(seq1);
      ProteinSequence protein2 = new ProteinSequence(seq2);
      //System.out.println(s1);
      //System.out.println(s2);
      // Step 2. do pairwised sequence alignement (Default is Local,i.e., Smith-Watermann).
      PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> seqAlignment =
        Alignments.getPairwiseAligner(protein1, protein2, params.alignType, params.gapPenalty, params.substitutionMatrix);

      // Step 3. find corresponding pairs.
      SequencePair<ProteinSequence, AminoAcidCompound> pair = seqAlignment.getPair();
      List<Atom> nCa1 = new ArrayList<Atom>();
      List<Atom> nCa2 = new ArrayList<Atom>();
      int nPairs = pair.getLength();
      int nIdentical = 0;
      int nSimilar   = 0;
      int nGaps      = 0;
      int numNonGapCols = 0;
      StringBuilder alnQuery  = new StringBuilder();
      StringBuilder alnTarget = new StringBuilder();
      StringBuilder alnSymbol = new StringBuilder();

      Compound gapSymbol =  AminoAcidCompoundSet.getAminoAcidCompoundSet().getCompoundForString("-");
      int[] align_se1 = new int[nPairs+1];
      int[] align_se2 = new int[nPairs+1];

      for (int i=1; i<=nPairs; i++){
        int currPos = i-1;
        Compound s1 = pair.getCompoundAt(1,i);
        Compound s2 = pair.getCompoundAt(2,i);
        //Position in original ca1 or ca2
        int Pos1 = pair.getIndexInQueryAt(i)-1;
        int Pos2 = pair.getIndexInTargetAt(i)-1;

        if ( !s1.equals(gapSymbol) && !s2.equals(gapSymbol) ){
          // normal amino acid pairs in the column
          numNonGapCols++;

          nCa1.add(ca1[Pos1]);
          nCa2.add(ca2[Pos2]);

          alnQuery.append(s1.getShortName());
          alnTarget.append(s2.getShortName());

          if(s1==s2){//An identical pair
            nIdentical++;
            nSimilar++;
            alnSymbol.append("|");
          }else{     //A similar pair
            nSimilar++;
            alnSymbol.append(":");
          }
          //Record the position in the array of ca1 or ca2
          align_se1[currPos] = Pos1;
          align_se2[currPos] = Pos2;

        }else if( s1.equals(gapSymbol) && s2.equals(gapSymbol) ){
          // There are all gaps in the column
          nGaps = nGaps+2;
          alnQuery.append("-");
          alnTarget.append(" ");
          alnSymbol.append(" ");
          align_se1[currPos] = -1;
          align_se2[currPos] = -1;
        }else {
          // There is a gap at least in the column
          nGaps++;
          alnSymbol.append(" ");

          align_se1[currPos] = -1;
          align_se2[currPos] = -1;

          if ( s1.equals(gapSymbol) ){
            alnQuery.append("-");
            alnTarget.append(s2.getShortName());

            align_se2[currPos] = Pos2;

          }else{
            alnQuery.append(s1.getShortName());
            alnTarget.append("-");

            align_se1[currPos] = Pos1;
          }
        }
      }
      // Print out the alignment
      //System.out.println(alnQuery.toString());
      //System.out.println(alnSymbol.toString());
      //System.out.println(alnTarget.toString());
      //printPSAinfo(seqAlignment);

      afpChain.setAlnLength(nPairs);
      afpChain.setGapLen(nGaps);
      afpChain.setAlnseq1(alnQuery.toString().toCharArray());
      afpChain.setAlnseq2(alnTarget.toString().toCharArray());
      afpChain.setAlnsymb(alnSymbol.toString().toCharArray());

      afpChain.setIdentity(nIdentical*1.0/numNonGapCols);
      afpChain.setSimilarity(nSimilar*1.0/numNonGapCols);
      afpChain.setOptLen(new int[]{numNonGapCols});
      afpChain.setAlignScore(seqAlignment.getScore());
      afpChain.setOptLength(numNonGapCols);
      // original length of ca1 and ca2
      afpChain.setCa1Length(ca1.length);
      afpChain.setCa2Length(ca2.length);
      //afpChain.setCa1Length(nCa1.size());
      //afpChain.setCa2Length(nCa2.size());

      //Less than 4 atoms
      if (alnQuery.length() < 4){
        return afpChain;
      }

      // do a 3D alignment...
      Atom[] qCa1 = nCa1.toArray(new Atom[0]);
      Atom[] tCa2 = nCa2.toArray(new Atom[0]);

      //System.out.println("Query  atoms: "+qCa1.length);
      //System.out.println("Target atoms: "+tCa2.length);

      // SVD to get Rotation coords and Shift coords
      SVDSuperimposer svds = new SVDSuperimposer(qCa1, tCa2);
      this.rotMatrix = svds.getRotation();
      this.tranMatrix = svds.getTranslation();

      // now we have all the info to perform the rotations ...
      for (Atom X : tCa2 ){
        Calc.rotate(X,this.rotMatrix);
        // shift Ca2 onto Ca1 ...
        Calc.shift(X,this.tranMatrix);
      }
      // Assert
      assert(qCa1.length == tCa2.length);
      double rmsd = SVDSuperimposer.getRMS(qCa1, tCa2);
      // Configure the  afpChain
      afpChain.setBlockRmsd(new double[]{rmsd});
      afpChain.setOptRmsd(new double[]{rmsd});
      afpChain.setTotalRmsdOpt(rmsd);
      afpChain.setChainRmsd(rmsd);

      //Using ce: refer to the sample SmithWaterman3Daligner.java
      CeParameters params = new CeParameters();
      CECalculator cecalc = new CECalculator(params);
      cecalc.setnAtom(numNonGapCols);
      //align_se1 and align_se2 are necessary for convertAfpChain(...)
      cecalc.setAlign_se1(align_se1);
      cecalc.setAlign_se2(align_se2);
      //cecalc.setMatMatrix(this.rotMatrix.getArray());
      cecalc.setLcmp(nPairs);
      // Create corresponding pairs from Structure info
      // utilize the routine in ce
      cecalc.convertAfpChain(afpChain, ca1, ca2);
      // Elasped time
      this.timeEnd = System.currentTimeMillis();
      afpChain.setIoTime(this.timeEnd-this.timeStart);
      // save the rotation matrix and the shift vector
      afpChain.setBlockRotationMatrix(new Matrix[]{this.rotMatrix});
      afpChain.setBlockShiftVector(new Atom[]{this.tranMatrix});


    } catch (Exception e) {
      throw new StructureException(e.getMessage(),e);
    }
    return afpChain;
  }

  @SuppressWarnings("unused")
private static char getOneLetter(Group g) {
    try {
      Character c = StructureTools.get1LetterCode(g.getPDBName());
      return c;
    } catch (Exception e){
      return 'X';
    }
  }

  @SuppressWarnings("unused")
private void printPSAinfo(PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> align){

    SequencePair<ProteinSequence, AminoAcidCompound> pair = align.getPair();
    int nPairs = pair.getLength();

    //Atom[] strBuf1 = new Atom[nPairs];
    //Atom[] strBuf2 = new Atom[nPairs];

    char[] alnseq1 = new char[nPairs];
    char[] alnseq2 = new char[nPairs] ;
    char[] alnsymb = new char[nPairs];

    //    int nIdentical = 0;
    //    int nSimilar   = 0;

    Compound gapSymbol =  AminoAcidCompoundSet.getAminoAcidCompoundSet().getCompoundForString("-");
    int sPos1 = pair.getIndexInQueryAt(1)-1;
    int ePos1 = pair.getIndexInQueryAt(nPairs)-1;

    int sPos2 = pair.getIndexInTargetAt(1)-1;
    int ePos2 = pair.getIndexInTargetAt(nPairs)-1;

    //System.out.print("Query: "+pos+" ");
    for (int i=1; i<=nPairs; i++){
      int currentPos = i-1;
      Compound s1 = pair.getCompoundAt(1,i);
      Compound s2 = pair.getCompoundAt(2,i);

      // normal amino acid pairs in the column
      alnseq1[i-1] = s1.getShortName().charAt(0);
      alnseq2[i-1] = s2.getShortName().charAt(0);
      if ( !s1.equals(gapSymbol) && !s2.equals(gapSymbol) ){
        if (s1 == s2){
          alnsymb[currentPos] = '|';
        }else{
          alnsymb[currentPos] = ':';
        }
      }else{
        alnsymb[currentPos] = ' ';
        //        if (s1.equals(gapSymbol)){

        //        }
      }
    }
    System.out.print("Query:  "+sPos1+" ");
    for (int i=1; i<=nPairs; i++){
      System.out.print(alnseq1[i-1]);
    }
    System.out.println(" "+ePos1);
    System.out.print("        ");
    for (int i=1; i<=nPairs; i++){
      System.out.print(alnsymb[i-1]);
    }
    System.out.println();
    System.out.print("Target: "+sPos2+" ");
    for (int i=1; i<=nPairs; i++){
      System.out.print(alnseq2[i-1]);
    }
    System.out.println(" "+ePos2);
    System.out.println("Max score: "+align.getMaxScore());
    System.out.println("Min score: "+align.getMinScore());
    System.out.println("Score: "+align.getScore());
    System.out.println("Similarity: "+align.getSimilarity());
  }

  public double[][] initSumOfDistances(int nse1, int nse2,Atom[] ca1, Atom[] ca2){
    double[][] mat = new double[nse1][nse2];
    return mat;
  }

}