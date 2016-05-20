Structure Comparison Tools for Java ver. 1.00
##############################################

:Authors: Wei Cao davecao@bi.a.u-tokyo.ac.jp  
:Version: 1.00 initialized on 5 Jan 2012  


I. Description
##################

Structure Comparison Tools for Java is developed for pairwise comparison of protein 3D structures, based on a open-source java framework for biological data (biojava). There are four comparison methods provided as following:

1. FATCAT rigid. (Ref.a)
2. FATCAT flexible. 
3. Combinatorial extenstion(Ce). (Ref.b)
4. CE circular permutation (CeCP), which sets default parameters to be appropriate for finding circular permutations.  

Reference:
a. Yuzhen Ye & Adam Godzik (2003) Flexible structure alignment by chaining aligned fragment pairs allowing twists. Bioinformatics vol.19 suppl. 2. ii246-ii255.   
  
b. Shindyalov IN, Bourne PE (1998) Protein structure alignment by incremental combinatorial extension (CE) of the optimal path. Protein Eng 11: 739-747 

.. ..
see also: (PDB comparison tools)
    Andreas Prlic; Spencer Bliven; Peter W. Rose; Wolfgang F. Bluhm; Chris Bizon; 
    Adam Godzik; Philip E. Bourne (2010)
    Pre-calculated protein structure alignments at the RCSB PDB website
    Bioinformatics 26: 2983-2985
  
II. Prerequisite:
##################

1. Apache Maven - A software project management and comprehension tool.
2. Java Runtime - version 1.6 or later 
3. Install the Jmol into local repository (If jmol is not available)

    prompt> cd path<--(root of source, see Source code tree)

    prompt> mvn install:install-file \
                -DgroupId=org.jmol \
                -DartifactId=jmol \
                -Dversion=12.0.22 \
                -Dpackaging=jar \
                -Dfile=./libs/Jmol.jar

III. Compilation:
##################

cd ROOT (ROOT is same location of pom.xml)  
    
mvn clean	 

mvn package  

The executable jar files can be located at 'jars' directory.
To run bilab-structure-1.0.jar, it needs the lib/*.jars on the java classpath.
On the contrary,  bilab-structure-1.0-jar-with-dependencies.jar can be run alone since all necessary libraries had been packaged into it.
i.e.,

    cmd> java -jar bilab-structure-1.0-jar-with-dependencies.jar [options]

IV. Usage:
##################

usage: java -jar bilab-structures-*.jar [options]

    --mol1, -a <Required>          The first molecule in PDB format.
    --alignAlgo <Optional>         The internally used pairwised alignment. Default is SW_local.
                                     **NW_global**: Needleman-Wunsch/Gotoh.
                                     **SW_local**: Smith-Waterman/Gotoh
                                     **GU_linear**: GUan_Uberbacher.
                                     **SW_linear**: Smith-Waterman/Gotoh with smart traceback
    --alignSeqRes                  If it presents, align the ATOM and
                                   SEQRES residues when parsing PDB files.
    --alignXMLfile <Optional>      Pairwised alignment file created from
                                   external program Blastp. If it is
                                   given, the internal pairwised sequence
                                   alignment will not be done.
    --mol2, -b <Required>          The second molecule in PDB format.

    --chain1 <Optional>         Specify the chain name of the first molecule in PDB format. 
                                   If not given, all alpha carbon atoms will be picked up 
                                   for pairwise structure alignment. (-c1 for short)
    --chain2 <Optional>         Specify the chain name of the second molecule in PDB format. Optioned.
                                   If not given, all alpha carbon atoms will be picked up
                                   for pairwise structure alignment.(-c2 for short)
    --showElapsedTime, -e          Print out elapsed time (boolean).
    --method, -m <Default=1>        Select comparison method(number):
                                        **1**. FATCAT rigid.
                                        **2**. FATCAT flexible.
                                        **3**. Combinatorial extenstion(CE).
                                        **4**. CE circular permutation(CECP).
                                        **5**. CE circular permutation side chain(CECPSideChain).
                                        **6**. Sequence-based comparison
    --showMemoryInfo, -mem             Print out used memory info(boolean).
    --output, -o<Optional>            The output file name.
    --parseCAonly                  If it presents, only CA atoms will be
                                     attained when parsing PDB files.
    --parseSecStruct               If it presents, parse secondary
                                     structures when parsing PDB files.
    --gapExt <Optional>    Gap Extension penalty for Sequence-based 
                                  structural alignment. (-ge) Default is 1
    --gapOpen <Optional>   Gap Open penalty for Sequence-based 
                                  structural alignment. (-go) Default is 5

    --outputFormat, -t <Default=xml>   The output file format:
                                     Raw format: raw.
                                     xml format: xml.
                                     nice summary: pretty.
    --gui, -g                          Show the pairwise comparison in graphic user interface.
    --using-gui, -u                  Do the pairwise comparison with a
                                     simple GUI. If this option is
                                     specified, others options will be
                                     ignored.
    --help, -h                         Print out usage.

e.g., chain A of 1CDG  v.s. chain B of 1TIM
 
    java -jar jars/bilab-structure-1.0-jar-with-dependencies.jar -a pdbs/1MI7.pdb -b pdbs/3WRP.pdb -c1 R -c2 A 


Result: the attributes in root node of the output xml 
::

method="jFatCat_rigid"  

probability="1.15e-01"   

alignScore="186.62"  

totalRmsdOpt="3.92"  

identity="0.0498"  

The above result is same as the pre-calculated results on the PDB site
http://www.rcsb.org/pdb/workbench/showPrecalcAlignment.do?action=pw_fatcat&name1=1CDG.A&name2=1TIM.B


IV. Run Jmol directly 
###########################

    java -classpath path/bilab-structure-1.0-jar-with-dependencies.jar org.openscience.jmol.app.Jmol


V. Run SimpleAlignmentGUI directly 
##################################

  java -jar jars/bilab-structure-1.0-jar-with-dependencies.jar -u

.. ..
Note for using the option --alignXMLfile with -m 6 (-m 6 means using the method, jSeqBase)

To use an external pre-existed pairwise-alignment file (blastp), 
you need to set -outfmt to use 5, i.e. produce results in the xml format. 
So far, this program can only read the xml output from blastp.

a. The program will use the first HSP segment to create the rotation matrix if there are several Hsp exists.

b. The program will terminate when it meets "No hit found" in the PSA alignment file generated by blastp.
 
