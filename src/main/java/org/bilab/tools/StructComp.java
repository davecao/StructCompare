package org.bilab.tools;

import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.*;
import org.biojava.bio.structure.align.*;
//import org.biojava.bio.structure.align.fatcat.*;
import org.biojava.bio.structure.align.fatcat.calc.*;
import org.biojava.bio.structure.align.model.*;
import org.biojava.bio.structure.align.ce.*;
import org.biojava.bio.structure.align.gui.*;
import org.biojava.bio.structure.align.xml.*;
//import org.biojava.bio.structure.align.util.*;
// apache commons cli
import org.apache.commons.cli.GnuParser;
//import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
//import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.OptionBuilder;
//import org.apache.commons.cli.OptionGroup;
//apache commons lang
import org.apache.commons.lang3.ArrayUtils;
// relection
import java.lang.reflect.*;

import java.io.*;
import java.util.concurrent.*;

import static java.lang.System.out;

// Bilab
import org.bilab.tools.SeqStructComp;
import org.bilab.tools.align.SscParameters;

public class StructComp
{
  private static final long MEGABYTE = 1024L * 1024L;

  private static String mol1 = null;//"1cdg.A";
  private static String chain1 = null;

  private static String mol2 = null;//"1tim.B";
  private static String chain2 = null;

  private static String cpMethod = "FatCatRigid";//"FR,FF,CE,CECP,CeSC"
  private static int cpMethodInx = 0;
  private static boolean showMemInfo = false;
  private static boolean showElapsedTime = false;
  private static boolean withGUI = false;
  private static boolean useGUI = false;
  private static boolean parseCAonly = false;
  private static boolean parseSecStruct = false;
  private static boolean alignSeqRes = false;
  // Parameters for Seqeuence-based structural alignment
  private static int gapOpen = 5;
  private static int gapExt  = 1;
  private static String blastXMLfile = "";
  private static String PSAalgo = "SW_local";
  private static String PSAalgorithms[] = {
    "NW_global",
    "SW_local",
    "GU_linear",
    "SW_linear"
  };

  private static String algorithms[] = {
    "FatCatRigid",
    "FatCatFlexible",
    "CeMain",
    "CeCPMain",
    "CeSideChainMain",
    "SeqStructComp"
  };

  private static String algoClassName[] = {
    "org.biojava.bio.structure.align.fatcat.FatCatRigid",
    "org.biojava.bio.structure.align.fatcat.FatCatFlexible",
    "org.biojava.bio.structure.align.ce.CeMain",
    "org.biojava.bio.structure.align.ce.CeCPMain",
    "org.biojava.bio.structure.align.ce.CeSideChainMain",
    "org.bilab.tools.SeqStructComp"
  };
  private static String outFormat[] = {"xml","raw","pretty"};
  private static Options opt = new Options();
  private static GnuParser parser = new GnuParser();
  private static CommandLine cl = null;
  private static HelpFormatter helpFormatter = new HelpFormatter();

  private static String oformat = "xml";
  //   private Boolean remoteFiles = false;
  //   private String pdbFilePath = "tmp";
  private static String ofName = null;

  public static long bytesToMegabytes(long bytes) {
    return bytes / MEGABYTE;
  }
  public static String printElaspedTime(long ms) {
    if (ms < 0) throw new IllegalArgumentException("Error: Time in milliseconds is negative.");
    StringBuilder sb = new StringBuilder(64);

    long days = TimeUnit.MILLISECONDS.toDays(ms);
    ms -= TimeUnit.MILLISECONDS.toMillis(days);

    long hours = TimeUnit.MILLISECONDS.toHours(ms);
    ms -= TimeUnit.MILLISECONDS.toMillis(hours);

    long minutes = TimeUnit.MILLISECONDS.toMinutes(ms);
    ms -= TimeUnit.MILLISECONDS.toMillis(minutes);

    long seconds = TimeUnit.MILLISECONDS.toSeconds(ms);

    sb.append(days);
    sb.append(" Days ");
    sb.append(hours);
    sb.append(" Hours ");
    sb.append(minutes);
    sb.append(" Minutes ");
    sb.append(seconds);
    sb.append(" Seconds ");

    return sb.toString();
  }

  // print help
  private static void printHelp() {
		//private HelpFormatter helpFormatter = new HelpFormatter();
    helpFormatter.printHelp("java -jar bilab-structures-*.jar",opt);
    System.exit(-1);
  }

  private static void printHelp(String message) {
    System.out.println(message);
    helpFormatter.printHelp("java -jar bilab-structures-*.jar ",opt);
    System.exit(-1);
  }

  // create options
  @SuppressWarnings("static-access")
  private void createOpts() {
    //Options opt = new Options();
    opt.addOption(
                  OptionBuilder.withLongOpt("method")
                  .withDescription("comparison method(number): \n"+
                                   "  1. FATCAT rigid.\n "+
                                   "  2. FATCAT flexible.\n "+
                                   "  3. Combinatorial extenstion(CE).\n "+
                                   "  4. CE circular permutation(CECP).\n" +
                                   "  5. CE circular permutation side chain(CECPSideChain).\n"+
                                   "  6. Sequence-based comparison")
                  .withType(Number.class)
                  .hasArg()
                  .withArgName("Default=1")
                  .create('m'));

    opt.addOption(
                  OptionBuilder.withLongOpt("mol1")
                  .withDescription("The first molecule in PDB format.")
                  .withType(String.class)
                  .hasArg()
                  .withArgName("Required")
                  //.isRequired(true)
                  .create('a'));

    opt.addOption(
                  OptionBuilder.withLongOpt("chain1")
                  .withDescription("Specify the chain name of the first molecule in PDB format.\n"+
                                   "If not given, all alpha carbon atoms will be picked up "+
                                   "for pairwise structure alignment.")
                  .withType(String.class)
                  .hasArg()
                  .withArgName("Optional")
                  .create("c1"));

    opt.addOption(
                  OptionBuilder.withLongOpt("mol2")
                  .withDescription("The second molecule in PDB format.")
                  .withType(String.class)
                  .hasArg()
                  .withArgName("Required")
                  //.isRequired(true)
                  .create('b'));

    opt.addOption(
                  OptionBuilder.withLongOpt("chain2")
                  .withDescription("Specify the chain name of the second molecule in PDB format. Optioned.\n"+
                                   "If not given, all alpha carbon atoms will be picked up "+
                                   "for pairwise structure alignment.")
                  .withType(String.class)
                  .hasArg()
                  .withArgName("Optional")
                  .create("c2"));
    opt.addOption(
                  OptionBuilder.withLongOpt("parseCAonly")
                  .withDescription("If it presents, only CA atoms will be attained when parsing PDB files.")
                  .hasArg(false)
                  .create());

    opt.addOption(
                  OptionBuilder.withLongOpt("parseSecStruct")
                  .withDescription("If it presents, parse secondary structures when parsing PDB files.")
                  .hasArg(false)
                  .create());

    opt.addOption(
                  OptionBuilder.withLongOpt("alignSeqRes")
                  .withDescription("If it presents, align the ATOM and SEQRES residues when parsing PDB files.")
                  .hasArg(false)
                  .create());

    opt.addOption(
                  OptionBuilder.withLongOpt("alignAlgo")
                  .withDescription("The internally used pairwised alignment.Default is SW_local.\n"
                                   +"NW_global: Needleman-Wunsch/Gotoh.\n"
                                   +"SW_local : Smith-Waterman/Gotoh\n"
                                   +"GU_linear: GUan_Uberbacher.\n"
                                   +"SW_linear: Smith-Waterman/Gotoh with smart traceback")
                  .withType(String.class)
                  .hasArg()
                  .withArgName("Optional")
                  .create());

    opt.addOption(
                  OptionBuilder.withLongOpt("gapOpen")
                  .withDescription("Gap Open penalty for Sequence-based structural alignment.Default is 5")
                  .withType(Integer.class)
                  .hasArg()
                  .withArgName("Optional")
                  .create("go"));

    opt.addOption(
                  OptionBuilder.withLongOpt("gapExt")
                  .withDescription("Gap Extension penalty for Sequence-based structural alignment.Default is 1")
                  .withType(Integer.class)
                  .hasArg()
                  .withArgName("Optional")
                  .create("ge"));

    opt.addOption(
                  OptionBuilder.withLongOpt("alignXMLfile")
                  .withDescription("Pairwised alignment file created from external program Blastp. If it is given, the internal pairwised sequence alignment will not be done.")
                  .withType(String.class)
                  .hasArg()
                  .withArgName("Optional")
                  .create());


    opt.addOption(
                  OptionBuilder.withLongOpt("output")
                  .withDescription("The output file name.")
                  .withType(String.class)
                  .hasArg()
                  .withArgName("Optional")
                  .create('o'));

    opt.addOption(
                  OptionBuilder.withLongOpt("outputFormat")
                  .withDescription("The output file format: \n"+
                                   "   Raw format: raw.\n"+
                                   "   xml format: xml.\n"+
                                   "   nice summary: pretty.")
                  .withType(String.class)
                  .hasArg()
                  .withArgName("Default=xml")
                  .create('t'));

    opt.addOption(
                  OptionBuilder.withLongOpt("showMemoryInfo")
                  .withDescription("Print out used memory info(boolean).")
                  //.withType(String.class)
                  .hasArg(false)
                  //.withArgName("Default=F")
                  .create("mem"));

    opt.addOption(
                  OptionBuilder.withLongOpt("showElapsedTime")
                  .withDescription("Print out elapsed time (boolean).")
                  //.withType(String.class)
                  .hasArg(false)
                  //.withArgName("Default=F")
                  .create("e"));
    opt.addOption(
                  OptionBuilder.withLongOpt("gui")
                  .withDescription("Show the pairwise comparison in graphic user interface.")
                  //.withType(String.class)
                  .hasArg(false)
                  //.withArgName("Default=F")
                  .create("g"));
    opt.addOption(
                  OptionBuilder.withLongOpt("using-gui")
                  .withDescription("Do the pairwise comparison with a simple GUI. If this option is specified, others options will be ignored.")
                  //.withType(String.class)
                  .hasArg(false)
                  //.withArgName("Default=F")
                  .create("u"));
    opt.addOption(
                  OptionBuilder.withLongOpt("help")
                  .withDescription("Print out usage.")
                  .hasArg(false)
                  .create("h"));
  }

  // Constructor
  private StructComp(String[] args) {
    // Create command line arguments
    createOpts();
    // Parse command line arguments
    try {
      cl = parser.parse(opt,args);
      if( cl.hasOption("using-gui") ){
        /*
        String ugui = cl.getOptionValue("using-gui");
        if( ugui.equals("T") ){
          useGUI = true;
        }else if( ugui.equals("F") ){
          //default has been set.
        }else{
          printHelp("Unrecognized the argument: --using-gui: "+ugui+".");
        }*/
        useGUI = true;
      }else{

        if( cl.hasOption("help") ){
          printHelp();
        }
        // Options: a and b
        if( !cl.hasOption("mol1") || !cl.hasOption("mol2")){
          //helpFormatter.printHelp("Usage",opt);
          printHelp();
        }else{
          mol1 = cl.getOptionValue("mol1");
          mol2 = cl.getOptionValue("mol2");
        }

        // Options: chain1 and chain2
        if( cl.hasOption("chain1") ){
          chain1 = cl.getOptionValue("chain1");
        }
        if( cl.hasOption("chain2") ){
          chain2 = cl.getOptionValue("chain2");
        }

        if( cl.hasOption("method") ){
          int i = ((Number)cl.getParsedOptionValue("method")).intValue() - 1;
          if( i>=0 && i<=5 ){
            cpMethod = algorithms[i];
            cpMethodInx = i;
          }else {
            printHelp("Unrecognized argument of --method: " + i+ ".");
          }
        }
        // The following 3 parameters for parsing PDB files
        if( cl.hasOption("parseCAonly") ){
          parseCAonly = true;
        }
        if( cl.hasOption("parseSecStruct") ){
          parseSecStruct = true;
        }
        if( cl.hasOption("alignSeqRes") ){
          alignSeqRes = true;
        }
        // Parameters for Sequence alignment
        if( cl.hasOption("alignAlgo") ){
          String temp = cl.getOptionValue("alignAlgo");
          boolean found = false;
          //Test whether it is one of NW_global, SW_local, GU_linear or SW_linear
          for( String x: PSAalgorithms ){
            if( x.contains(temp) ) {
              found = true;
              PSAalgo = temp;
            }
          }
          if (!found) {
        	  throw new org.apache.commons.cli.UnrecognizedOptionException("The PSA algorithm is unknown."+temp);
          }
        }

        if( cl.hasOption("gapOpen") ){
          gapOpen = ((Number)cl.getParsedOptionValue("gapOpen")).intValue();
        }
        if( cl.hasOption("gapExt") ){
          gapExt = ((Number)cl.getParsedOptionValue("gapExt")).intValue();
        }
        if( cl.hasOption("alignXMLfile") ){
          String temp = cl.getOptionValue("alignXMLfile");
          //Test the existence of the alignment file.
          boolean exists = (new File(temp)).exists();
          if (exists) {
              blastXMLfile = temp;
          } else {
        	  System.out.println("The Specified option: --alignXMLfile.");
        	  System.out.println("Could not find the alignment file: "+temp);
        	  throw new java.io.FileNotFoundException("Could not find the alignment file:"+temp);
          }
        }
        if( cl.hasOption("showMemoryInfo") ){
          /*
          String si = cl.getOptionValue("showMemoryInfo");
          if( si.equals("T") ){
            showMemInfo = true;
          }else if( si.equals("F") ){
            //default has been set.
          }else{
            printHelp("Unrecognized the argument: --showMemoryInfo: "+si+".");
          }
          */
          showMemInfo = true;
        }
        if( cl.hasOption("showElapsedTime") ){
          /*
          String st = cl.getOptionValue("showElapsedTime");
          if( st.equals("T") ){
            showElapsedTime = true;
          }else if( st.equals("F") ){
            //default has been set.
          }else{
            printHelp("Unrecognized the argument of --showElapsed: "+st+".");
          }
          */
          showElapsedTime = true;
        }
        /*
          if( this.cl.hasOption("r") && this.cl.hasOption("p") ){
          this.remoteFiles = this.cl.getOptionValue("r");
          this.pdbFilePath = this.cl.getOptionValue("p");
          }
        */
        if( cl.hasOption("output") ){
          ofName = cl.getOptionValue("output");
        }

        if( cl.hasOption("gui") ){
          /*
          String gui = cl.getOptionValue("gui");
          if( gui.equals("T") ){
            withGUI = true;
          }else if( gui.equals("F") ){
            //default has been set.
          }else{
            printHelp("Unrecognized the argument: --gui: "+gui+".");
          }
          */
          withGUI = true;
        }

        if( cl.hasOption("outputFormat") ){
          oformat = cl.getOptionValue("outputFormat");
          int indexOfformat = ArrayUtils.indexOf(outFormat,oformat);
          if( indexOfformat == -1 ) {
            printHelp("Unsupport output format"+oformat);
          }
        }
      }

    }catch(org.apache.commons.cli.ParseException e) {
      //helpFormatter.printHelp("Usage",opt);
      printHelp("Error CMD arguments:"+e.getMessage());
    } catch (FileNotFoundException e) {
		// TODO Auto-generated catch block
		//e.printStackTrace();
		System.exit(0);
	}
  }
  /**
   * Static factory method: it can be called directly.
   * @param command line arguments
   */
  public static StructComp newInstance(String[] args){
    return new StructComp(args);
  }

  /**
   * Get CA atoms from Structure s by the specified chain Id
   * @Description extract CA atoms from the structure object
   * @param s
   * @param ch
   */
  private static Atom[] getAlignedCAAtoms(Structure s, String ch) throws org.biojava.bio.structure.StructureException {
    Atom[] ca = null;
    Chain chain = null;
    java.util.ArrayList<Atom> caList =new java.util.ArrayList<Atom>();
    if( ch == null || ch.equals("") ){// if the chain name is not specified
      ca = StructureTools.getAtomCAArray(s); // include HETATM either
    }else{
      if( s.hasChain(ch) ){// Exist of the specified chain?
        chain = s.findChain(ch);
        // Chain may include hetatm either.
        //ca = StructureTools.getAtomCAArray(chain);
        java.util.List<Group> groups = chain.getAtomGroups();

        for( Group g: groups ){
          //Group type: amino, hetatm and nucleotide.
          String gtp = g.getType();
          if( gtp.equals("amino")){
            if(g.hasAtom("CA")){
                caList.add(g.getAtom("CA"));
            }else{
               System.err.println(s.getPDBCode()+":Residue "+g.getPDBName()+" "+g.getResidueNumber().getSeqNum()+" of "+g.getChainId());
            }
          }
        }
        ca = caList.toArray(new Atom[]{});
      }else{
        System.out.println("Could not find the specified chain "+ch+" in "+s.getPDBCode());
        System.exit(-1);
      }
    }
    //out.format("Group type: %s\n",gtp);
    //out.format("%s %d\n",s.getPDBCode(),ca.length);
    /*
      for (Atom x: ca){
      out.format("%d %s %s %s\n",x.getPDBserial(),x.getFullName(),x.getGroup().getChainId(),x.getGroup().getPDBName());
      }
    */
    return ca;
  }

  private static String getAlgoNameFromString (String algoName){
    String s = null;
    String fieldname = "algorithmName";
    Class<?> algo;
    //    Object obj = null;
    try {
      algo = Class.forName(algoName);
      //obj = algo.newInstance();
      Field a = algo.getDeclaredField(fieldname);
      s = (String) a.get(fieldname);
      //System.out.println("Algorithm Name: "+ s);

    }catch(NoSuchFieldException e){
      System.out.println(e);
    }catch(IllegalAccessException e){
      System.out.println(e);
    }catch(SecurityException e){
      System.out.println(e);
    }catch(ClassNotFoundException e){
      throw new RuntimeException(e);
    }
    return s;
  }

  // get default parameters
  private static Object getAlgoDefaultParms(String algo) {
    if ( algo.equals("FatCatRigid")
         || algo.equals("FatCatFlexible")){
      return new FatCatParameters();
    }else if ( algo.equals("CeMain")
               || algo.equals("CeCPMain")
               ||  algo.equals("CeSideChain") ){
      return new CeParameters();
    }else if ( algo.equals("SeqStructComp")){
      SscParameters p =  new SscParameters();
      p.setGapOpen((short) gapOpen);
      p.setGapExtension((short) gapExt);
      p.setBlastAlignmentFile(blastXMLfile.split("#"));
      p.setAlgorithm(PSAalgo.split("#"));
      return p;
    }
    return null;
  }

  private static void showResult(AFPChain afpChain, Atom[] ca1, Atom[] ca2, String algo,String format) {
    try{
      if( format.equals("raw") ) {
        if ( algo.equals("FatCatRigid")
             || algo.equals("FatCatFlexible")){
          // show original FATCAT output:
          System.out.println(afpChain.toFatcat(ca1,ca2));
        }else if ( algo.equals("CeMain")
                   || algo.equals("CeCPMain")
                   || algo.equals("CeSideChain") 
                   || algo.equals("SeqStructComp")){
          // show original CE output:
          System.out.println(afpChain.toCE(ca1,ca2));
        }
      } else if ( format.equals("xml") ){
        // print XML representation
        System.out.println(AFPChainXMLConverter.toXML(afpChain,ca1,ca2));

      } else if ( format.equals("pretty") ){
        // show a nice summary print
        System.out.println(AfpChainWriter.toWebSiteDisplay(afpChain, ca1, ca2));
      }
      // print rotation matrices

      //System.out.println(afpChain.toRotMat());
      //System.out.println(afpChain.toCE(ca1, ca2));

      //StructureAlignmentDisplay.display(afpChain, ca1, ca2);
    }catch(Exception e) {
      e.printStackTrace();
      return;
    }
  }
  // output to file
  private static void showResult(String ofname,AFPChain afpChain,
                                 Atom[] ca1, Atom[] ca2,
                                 String algo,String format) {
    BufferedWriter o = null;

    try{
      File file = new File(ofname);
      boolean exist = file.createNewFile();
      if(!exist) {
        System.out.println(ofname+" will be overwrited.");
      }

      FileWriter fstream = new FileWriter(file,false);
      o = new BufferedWriter(fstream);

      if( format.equals("raw") ) {
        if ( algo.equals("FatCatRigid")
             || algo.equals("FatCatFlexible")){
          // show original FATCAT output:
          o.write(afpChain.toFatcat(ca1,ca2));

        }else if ( algo.equals("CeMain")
                   || algo.equals("CeCPMain")
                   ||  algo.equals("CeSideChain") ){
          // show original CE output:
          o.write(afpChain.toCE(ca1,ca2));
        }
      } else if ( format.equals("xml") ){
        // print XML representation
        o.write(AFPChainXMLConverter.toXML(afpChain,ca1,ca2));

      } else if ( format.equals("pretty") ){
        // show a nice summary print
        o.write(AfpChainWriter.toWebSiteDisplay(afpChain, ca1, ca2));
      }
    }catch(IOException e){
      System.err.println("Caught IOException: " + e.getMessage());
    }catch(Exception e) {
      e.printStackTrace();
      System.exit(-1);

    }finally {
      try {
        if (o != null){
          System.out.println("Write to file finished.");
          o.close();
        }else{
          System.out.println("FileWriter not open");
        }
      }catch(IOException e){
        System.err.println("Caught IOException: " + e.getMessage());
      }
    }
  }

  private static void printInfo(Structure s1, Structure s2) {
    // print info
    String pdbCode1 = s1.getPDBCode();
    String pdbCode2 = s2.getPDBCode();

    if(( chain1 != null ) && ( chain2 != null )){
      out.format("Method: %s\n Mol1:%s.%s Mol2:%s.%s\n output format:%s\n",
                 cpMethod,pdbCode1,chain1,pdbCode2,chain2,oformat);
    }else{
      out.format("Method: %s\n Mol1:%s Mol2:%s\n output format:%s\n",
                 cpMethod,pdbCode1,pdbCode2,oformat);
    }
  }

  public static void main( String[] args ){
    // for calculating elapsed tiem
    long startTime = 0L;

    // Using GUI interface
    //System.setProperty("PDB_DIR","/tmp/");
    //AlignmentGui.getInstance();

    //String pdbFilePath=args[0];
    //StructComp comparer = new StructComp(args);

    //Create an object by the static factory method in Class StructComp
    StructComp.newInstance(args);

    if ( showElapsedTime ) {
      startTime = System.currentTimeMillis();
    }
    //boolean isSplit = true;
    // Download from web
    //AtomCache cache = new AtomCache(this.pdbFilePath, isSplit);

    PDBFileReader pdbreader = new PDBFileReader();

    // configure the parameters of file parsing
    FileParsingParameters params = new FileParsingParameters();
    // Only parse CA
    if (parseCAonly){
      params.setParseCAOnly(true);
    }
    // should secondary structure get parsed from the file
    if (parseSecStruct) {
      params.setParseSecStruc(false);
    }
    // should the ATOM and SEQRES residues be aligned when creating the internal data model?
    if (alignSeqRes) {
      params.setAlignSeqRes(false);
    }
    //out.format("Max atoms:%d\n",params.getMaxAtoms());
    // set parser's parameter
    pdbreader.setFileParsingParameters(params);

    //pdbreader.setAlignSeqRes(false);
    Structure structure1 = null;
    Structure structure2 = null;

    try {
      StructureAlignmentFactory.addAlgorithm(new SeqStructComp());
      StructureAlignmentFactory.removeAlgorithm("Smith-Waterman superposition");
      if ( useGUI ){
        AlignmentGui.getInstance();
      }else{

        // To run FATCAT in the flexible variant say
        // FatCatFlexible.algorithmName below
               
        String algoName = getAlgoNameFromString(algoClassName[cpMethodInx]);
        StructureAlignment algorithm  = StructureAlignmentFactory.getAlgorithm(algoName);

        // use web
        //structure1 = cache.getStructure(this.mol1);
        //structure2 = cache.getStructure(this.mol2);

        // read pdb locally
        structure1 = pdbreader.getStructure(mol1);
        structure2 = pdbreader.getStructure(mol2);

        printInfo(structure1, structure2);
        // If options for chains are not specified,
        // all alpha carbon atoms will be picked up.
        //ca1 = StructureTools.getAtomCAArray(structure1);
        //ca2 = StructureTools.getAtomCAArray(structure2);

        Atom[] ca1 = getAlignedCAAtoms(structure1, chain1);
        Atom[] ca2 = getAlignedCAAtoms(structure2, chain2);

        // get default parameters
        //FatCatParameters params = new FatCatParameters();
        out.format("Do structure alginment... please wait\n");
        AFPChain afpChain = algorithm.align(ca1,ca2,getAlgoDefaultParms(cpMethod));
        afpChain.setName1(mol1);
        afpChain.setName2(mol2);

        //show Result
        if( ofName != null ){
          //output to file
          out.format("Trying to save to %s \n",ofName);
          showResult(ofName,afpChain, ca1, ca2, cpMethod,oformat);
        }else{
          // print to screen
          showResult(afpChain, ca1, ca2, cpMethod,oformat);
        }
        out.format("Calculation FINISHED.\n");

        if( showMemInfo ){
          // Get java runtime
          Runtime runtime = Runtime.getRuntime();
          // Run tje garbage collector
          runtime.gc();
          //calculate the used memory
          long memory = runtime.totalMemory() - runtime.freeMemory();
          out.format("Used memory: %s M [%s bytes]\n",bytesToMegabytes(memory),memory);
        }

        if( showElapsedTime ) {
          long stopTime = System.currentTimeMillis();
          long elapsedTime = stopTime-startTime;
          //          String e = null;

          out.format("The elapsed time: %s\n",printElaspedTime(elapsedTime));
        }
        if( withGUI ) {
          if ( afpChain.getBlockRotationMatrix() == null || afpChain.getBlockRotationMatrix().length == 0) {
            // probably the alignment is too short!
            System.err.println("No rotation matrix found!");
            System.err.println("The alignment may be too short!");
          }else{
              StructureAlignmentDisplay.display(afpChain,ca1,ca2);
          }
        }
      }
    } catch (Exception e) {
      e.printStackTrace();
      return;
    }
  }//end of main
}
