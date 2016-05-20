package org.bilab.tools.align;


import java.lang.reflect.Field;
//import java.lang.reflect.ParameterizedType;
//import java.lang.reflect.Type;
//import java.lang.reflect.TypeVariable;

import java.util.ArrayList;
//import java.util.Collections;
//import java.util.LinkedHashMap;
import java.util.List;
import java.util.logging.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import org.biojava3.core.util.XMLHelper;

import org.biojava3.alignment.SimpleAlignedSequence;
import org.biojava3.alignment.SimpleSequencePair;
import org.biojava3.alignment.template.AbstractPairwiseSequenceAligner;
import org.biojava3.alignment.template.AlignedSequence;
import org.biojava3.alignment.template.AlignedSequence.Step;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.template.SequencePair;
//import org.biojava3.alignment.template.PairwiseSequenceAligner;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
//import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;
//import org.biojava3.core.sequence.template.Sequence;

//@SuppressWarnings("rawtypes")
//public class BlastXMLParser<S extends Sequence<C>,C extends Compound> 
//	extends AbstractPairwiseSequenceAligner<S,C>
//	implements PairwiseSequenceAligner<S,C>  {
public class BlastXMLParser extends AbstractPairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> {
	private static final Logger log = Logger.getLogger(BlastXMLParser.class.getName());
	
	protected Document blastDoc = null;
//	private Class<S> sClass;
//	private Class<C> cClass;
	
//	private S targetSequence;
//	private S querySequence;
//	private C compund;
	private ProteinSequence orgTargetSequence;
	private ProteinSequence orgQuerySequence;
	
	private ProteinSequence targetSequence;
	private ProteinSequence querySequence;
//	private AminoAcidCompound compound;
	
	//Blast parameters
	protected String ScoringMatrix;
	protected double expect;
	protected double gapOpen;
	protected double gapExtend;
	protected boolean useFilter;
	
	// Query info
	protected String queryDef;
    protected String querySequenceStr;
	protected int queryLen;
	protected int queryFrom;
	protected int queryTo;

	// Hit info
	protected String hitId;
	protected String hitDef;
	protected String targetSequenceStr;
	protected int hitLen;
	
	protected int hitFrom;
	protected int hitTo;
	
	//Blast Hsp results
	protected double hspScore;
	protected double hspEvalue;
	protected double hspBitScore;
	protected int hspIdentity;
	protected int hspPositive;
	protected int hspGaps;
	protected int hspAlignLen;
	// No Hit found ?
	protected boolean NoHitFound = false;
	
	public static String xmlFile ="";
	
	
	/*
	public <S extends Sequence<C>,C extends Compound> BlastXMLParser<S,C> parse(){
	
		return parseBlastXML(xmlFile);
	}
	*/
	
	private BlastXMLParser(String blastFileName,ProteinSequence org1,ProteinSequence org2){
		xmlFile = blastFileName;
		orgQuerySequence = org1;
		orgTargetSequence = org2;
		
		log.info("Reading "+blastFileName);
		try {
			
			this.blastDoc = XMLHelper.loadXML(blastFileName);
			getBlastOutputParam();
			getBlastQueryInfo();
			getBlastHitInfo();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		log.info("Finished.");
	}
	
	//Static factory method for a direct call 
	//@SuppressWarnings("rawtypes")
	/*public static <S extends Sequence<C>,C extends Compound> BlastXMLParser<S,C> parseBlastXML(String blastFileName) {
		return new BlastXMLParser<S,C>(blastFileName);
	}
	*/
	public static BlastXMLParser parseBlastXML(String blastFileName,ProteinSequence org1,ProteinSequence org2) {
		return new BlastXMLParser(blastFileName,org1,org2);
	}
	/**
	 * Parsing parameters of blastp search
	 * 
	 * @throws Exception
	 */
	public void getBlastOutputParam() throws Exception{
		ArrayList<Element> elementList = XMLHelper.selectElements(
				this.blastDoc.getDocumentElement(),
				"BlastOutput_param/Parameters");
		for (Element element : elementList){
			Element pMatrix  = XMLHelper.selectSingleElement(element, "Parameters_matrix");
			Element pExpect  = XMLHelper.selectSingleElement(element, "Parameters_expect");
			Element pGapOpen = XMLHelper.selectSingleElement(element, "Parameters_gap-open");
			Element pGapExt  = XMLHelper.selectSingleElement(element, "Parameters_gap-extend");
			Element pFilter  = XMLHelper.selectSingleElement(element, "Parameters_filter");
			this.ScoringMatrix = pMatrix.getTextContent();
			this.expect        = Double.parseDouble(pExpect.getTextContent());
			this.gapOpen       = Double.parseDouble(pGapOpen.getTextContent());
			this.gapExtend     = Double.parseDouble(pGapExt.getTextContent());
			if (pFilter.getTextContent().equals("F")) {
				useFilter = false;
			}else{
				useFilter = true;
			}
		}
	}
	
	/**
	 *  Parsing the query length and the query description
	 *  
	 * @throws Exception
	 */
	public void getBlastQueryInfo() throws Exception{
		ArrayList<Element> elementList = XMLHelper.selectElements(
				this.blastDoc.getDocumentElement(),
				"BlastOutput_iterations/Iteration");
		for (Element element : elementList){
			//Element qId  = XMLHelper.selectSingleElement(element, "Iteration_query-ID");
			Element qDef = XMLHelper.selectSingleElement(element, "Iteration_query-def");
			Element qLen = XMLHelper.selectSingleElement(element, "Iteration_query-len");
			this.queryDef = qDef.getTextContent();
			//String tempLen = qLen.getTextContent();
			//this.queryLen = Integer.parseInt(tempLen);
			this.queryLen = Integer.parseInt(qLen.getTextContent());
		}
	}
	
	/**
	 * Get Hit info
	 * 
	 * @throws Exception
	 */
	public void getBlastHitInfo() throws Exception {
		//log.info("Query for Iteration hits");
		ArrayList<Element> IterationHitsList = XMLHelper.selectElements(
				this.blastDoc.getDocumentElement(),
				"BlastOutput_iterations/Iteration/Iteration_hits");
		log.info(IterationHitsList.size()+" iteration hits");
		if (IterationHitsList.size() == 0) {
		   this.NoHitFound();	
		}else{
			// get the first Iteration_hits
			Element firstIterationHitsElement = IterationHitsList.get(0);
			// find the first Hit in the Iteration_hits
			ArrayList<Element> hitList = XMLHelper.selectElements(firstIterationHitsElement, "Hit");
			if(hitList.size() == 0) {
				this.NoHitFound();
			}else{
				// get the first Hit
				Element firstHitElement = hitList.get(0);
				Element hHitId  = XMLHelper.selectSingleElement(firstHitElement, "Hit_id");
				Element hHitLen = XMLHelper.selectSingleElement(firstHitElement, "Hit_len");
				Element hHitDef = XMLHelper.selectSingleElement(firstHitElement, "Hit_def");
				this.hitId = hHitId.getTextContent();
				this.hitLen = Integer.parseInt(hHitLen.getTextContent());
				this.hitDef = hHitDef.getTextContent();
				// Process Hsps in the first Hit
				// get the first Hsp segment
				ArrayList<Element> hspList = XMLHelper.selectElements(
						firstHitElement, "Hit_hsps/Hsp");
				if(hspList.size() ==0){
					this.NoHitFound();
				}else{
					Element hsp = hspList.get(0);
					Element hHitFrom = XMLHelper.selectSingleElement(hsp, "Hsp_hit-from");
					Element hHitTo   = XMLHelper.selectSingleElement(hsp, "Hsp_hit-to");
					Element hQueryFrom = XMLHelper.selectSingleElement(hsp, "Hsp_query-from");
					Element hQueryTo   = XMLHelper.selectSingleElement(hsp, "Hsp_query-to");
					this.hitFrom = Integer.parseInt(hHitFrom.getTextContent());
					this.hitTo   = Integer.parseInt(hHitTo.getTextContent());
					this.queryFrom = Integer.parseInt(hQueryFrom.getTextContent());
					this.queryTo   = Integer.parseInt(hQueryTo.getTextContent());
					// Hsp scores
					Element pBitScore = XMLHelper.selectSingleElement(hsp, "Hsp_bit-score");
					Element pScore    = XMLHelper.selectSingleElement(hsp, "Hsp_score");
					Element pEvalue   = XMLHelper.selectSingleElement(hsp, "Hsp_evalue");
					Element pIdentity = XMLHelper.selectSingleElement(hsp, "Hsp_identity");
					Element pPositive = XMLHelper.selectSingleElement(hsp, "Hsp_positive");
					Element pGaps     = XMLHelper.selectSingleElement(hsp, "Hsp_gaps");
					Element pAlignLen = XMLHelper.selectSingleElement(hsp, "Hsp_align-len");
					Element pQuerySeq = XMLHelper.selectSingleElement(hsp, "Hsp_qseq");
					Element pTargetSeq= XMLHelper.selectSingleElement(hsp, "Hsp_hseq");
					
					this.hspBitScore = Double.parseDouble(pBitScore.getTextContent());
					this.hspScore    = Double.parseDouble(pScore.getTextContent());
					this.hspEvalue   = Double.parseDouble(pEvalue.getTextContent());
					this.hspIdentity = Integer.parseInt(pIdentity.getTextContent());
					this.hspPositive = Integer.parseInt(pPositive.getTextContent());
					this.hspGaps     = Integer.parseInt(pGaps.getTextContent());
					this.hspAlignLen = Integer.parseInt(pAlignLen.getTextContent());

					// Store sequence (Need java v.1.7)
					//this.querySequence = sClass.getConstructor(String.class).newInstance(pQuerySeq.getTextContent());
					//this.targetSequence = sClass.getConstructor().newInstance(pTargetSeq.getTextContent());
					//this.querySequence = this.setInstance(sClass, pQuerySeq.getTextContent());

					//getSInstance(sClass, querySequenceStr);
					//this.sClass = getSTypeParameterClass();
					//this.cClass = getCTypeParameterClass();
					//this.querySequence = getSTypeParameterClass().newInstance();
					this.querySequenceStr  = pQuerySeq.getTextContent();
					this.targetSequenceStr = pTargetSeq.getTextContent();
					this.querySequence = new ProteinSequence(querySequenceStr);
					this.targetSequence = new ProteinSequence(targetSequenceStr);
					//this.querySequence.setDescription(this.queryDef);
					//this.targetSequence.setDescription(this.hitDef);
				}
			}
			
		}
	}
	
	private void NoHitFound() {
		// No Hit found is true
		this.NoHitFound =true;
		// Could not find the corresponding pairs so as not to create
		// the rotation matrix and the shift vector
		System.out.println("No hit found from the input alignment file,"+xmlFile);
		System.out.println("Could not found the corresponding pairs so as not to create the rotation matrix.");
		System.exit(0);
	}
/*	public Class<S> getTypeParameterClass(){
		Type type = getClass().getGenericSuperclass();
		ParameterizedType paramType = (ParameterizedType) type;
		
	}
*/
	
/*	
	@SuppressWarnings("unchecked")
	public Class<S> getSTypeParameterClass(){
		Type type = getClass().getGenericSuperclass();
		if (type instanceof ParameterizedType){
			//ParameterizedType paramType = (ParameterizedType) ((Class)((ParameterizedType) type).getRawType()).getGenericSuperclass();
			ParameterizedType paramType = (ParameterizedType)type;
			Type tp0 = paramType.getActualTypeArguments()[0];
			return (Class<S>)tp0;
		}else if (type instanceof Class){
			return (Class<S>) type;
		}
		return null;
	}
	
	@SuppressWarnings("unchecked")
	public Class<C> getCTypeParameterClass(){
		Type type = this.getClass().getGenericSuperclass();
		if (type instanceof ParameterizedType){
			ParameterizedType paramType = (ParameterizedType) type;
			return (Class<C>)paramType.getActualTypeArguments()[1];
		}else if ( type instanceof Class){
			return (Class<C>)type;
		}
		return null;
	}
	public long getComputationTime() {
		// TODO Auto-generated method stub
		return 0;
	}
*/
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public Profile getProfile() {
		// TODO Auto-generated method stub
		return null;
	}

	public double getDistance() {
		// TODO Auto-generated method stub
		return 0;
	}

	public double getDistance(double arg0) {
		// TODO Auto-generated method stub
		return 0;
	}

	public int getMaxScore() {
		// TODO Auto-generated method stub
		return 0;
	}

	public int getMinScore() {
		// TODO Auto-generated method stub
		return 0;
	}

	public int getScore() {
		// TODO Auto-generated method stub
		return (int) this.hspScore;
	}

	public double getSimilarity() {
		// TODO Auto-generated method stub
		return 0;
	}

	public double getSimilarity(double arg0) {
		// TODO Auto-generated method stub
		return 0;
	}

	public ProteinSequence getQuery() {
		// TODO Auto-generated method stub
		return this.querySequence;
	}

	public ProteinSequence getTarget() {
		// TODO Auto-generated method stub
		return this.targetSequence;
	}
/*
	public SequencePair<S, C> getPair() {
		// TODO Auto-generated method stub
		return null;
	}
*/	
	@SuppressWarnings("unused")
	private int getCapacity(List<Step> qlist) throws Exception{
		Field dataField = ArrayList.class.getDeclaredField("elementData");
		dataField.setAccessible(true);
		return ((Object[])dataField.get(qlist)).length;
	}
	
	public SequencePair<ProteinSequence,AminoAcidCompound> getPair(){
		// TODO Auto-generated method stub
		int nPairs = this.querySequence.getLength();
		
		//List<AlignedSequence.Step> qlist = Collections.synchronizedList(new ArrayList<AlignedSequence.Step>(nPairs));
		//List<AlignedSequence.Step> tlist = Collections.synchronizedList(new ArrayList<AlignedSequence.Step>(nPairs));
		List<AlignedSequence.Step> qlist = new ArrayList<AlignedSequence.Step>(nPairs);
		List<AlignedSequence.Step> tlist = new ArrayList<AlignedSequence.Step>(nPairs);
		/*try {
			System.out.println("initial capacity of qlist:"+ getCapacity(qlist));
			System.out.println("initial capacity of tlist:"+ getCapacity(tlist));
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		*/
		Compound gapSymbol =  AminoAcidCompoundSet.getAminoAcidCompoundSet().getCompoundForString("-");
		
		if (this.querySequenceStr.length() != this.targetSequenceStr.length()){
			System.out.println("The length of aligned two segment are not same.");
			System.out.println("Check the alignment file: "+ xmlFile);
			System.exit(-1);
		}else{
			
			for (int i=0; i<nPairs; i++){
        		
				qlist.add(i, AlignedSequence.Step.GAP);
        		tlist.add(i, AlignedSequence.Step.GAP);
        		
				Compound s1 = querySequence.getCompoundAt(i+1);
		        Compound s2 = targetSequence.getCompoundAt(i+1);
		        
		        if ( !s1.equals(gapSymbol) && !s2.equals(gapSymbol) ){
		            // An identical or similar position 
		        	qlist.set(i, AlignedSequence.Step.COMPOUND);
		        	tlist.set(i, AlignedSequence.Step.COMPOUND);
		        }else{
		        	// At least has one gap
		        	if( ! s1.equals(gapSymbol) ){
		        		qlist.set(i, AlignedSequence.Step.COMPOUND);
		            }
		        	if( ! s2.equals(gapSymbol) ){
		        		tlist.set(i, AlignedSequence.Step.COMPOUND);
		            }
		        }
			}
		}
		//System.out.println("length of qlist: "+qlist.size());
		AlignedSequence<ProteinSequence,AminoAcidCompound> qAligned 
		   = new SimpleAlignedSequence<ProteinSequence, AminoAcidCompound>
		       (this.orgQuerySequence, qlist, this.queryFrom-1,this.orgQuerySequence.getLength()-this.queryTo);
		AlignedSequence<ProteinSequence,AminoAcidCompound> tAligned 
		   = new SimpleAlignedSequence<ProteinSequence, AminoAcidCompound>
		       (this.orgTargetSequence, tlist,this.hitFrom-1,this.orgTargetSequence.getLength()-this.hitTo);
		
		SequencePair<ProteinSequence,AminoAcidCompound> pair 
		   = new SimpleSequencePair<ProteinSequence, AminoAcidCompound>
		            (qAligned,tAligned);
		
		return pair;
	}
	
	@Override
	protected void setProfile(List<Step> arg0, List<Step> arg1) {
		// TODO Auto-generated method stub
		
	}
	
	/*
	 * @return E value from the pre-existed alignment file by Blastp  
	 */
	private double setToparamObj(){
		
		return this.hspEvalue;
	}
}