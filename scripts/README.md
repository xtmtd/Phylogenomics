
# Scripts of phologenomics

## Introduction




All ten bash scripts are tested in Linux systems. Their functions are described below.

1. BUSCO_extraction.sh: this script extracts the amino acid and nucleotides sequences for all the BUSCO results (run_* folders) by using TransDecoder, i.e., modifying the head name of the fasta files for each locus and merging sequences of the same locus into the fasta files, then filtering loci having too few taxa (less than three).
2. align_MAFFT.sh: this script performs multiple sequence alignment by using MAFFT or MAGUS.
3. trmming_alignments.sh: this script performs trimming of multiple sequence alignments by using trimAl, BMGE or ClipKIT, i.e., identifying and removing highly divergent sites, or retaining parsimony-informative/constant sites from multiple sequence alignments.
4. loci_filtering_alignment-based.sh: this script filters loci by detecting alignment length, number of parsimony-informative sites, percentage of parsimony-informative sites in the alignment, GC content, compositional heterogeneity (RCV, Relative composition variability), evolutionary rate (average pairwise identity) and likelihood mapping criteria by using PhyKIT or IQ-TREE.
5. gene_trees.sh: this script generates individual gene trees for a set of alignments by using IQ-TREE.
6. loci_filtering_tree-based.sh: this script filter loci using gene tree-based methods, i.e., average bootstrap support (ABS), Degree of violation of the molecular clock (DVMC), treeness, signal-to-noise ratio (treeness over rcv), and spurious homologs (possible paralogs, incorrectly assembled sequences) by using PhyKIT or IQ-TREE.
7. treeshrink.sh: this script performs detection of outlier long branches in collections of phylogenetic trees.
8. matrix_generation.sh: this script generates the supermaxtrix, partition and occupancy for alignments by using PhyKIT.
9. astral.sh: this script generates the species tree by using ASTRAL.
10. mcmctree_AA.sh: this script estimates divergence time by using MCMCTree for large amino acid dataset.

## Requirements

Some bioinformatic tools are neccessary for above scripts. Most of them are recommended to be added into the environmental paths. Softwares, versions and source are listed below. All the softwares must be installed in Linux systems.

   FASconCAT-G v1.04 (https://github.com/PatrickKueck/FASconCAT-G)  
   PHYLUCE v1.6.6 (http://phyluce.readthedocs.io/en/latest/index.html)  
   PhyKIT (https://github.com/JLSteenwyk/PhyKIT)  
   TransDecoder v5.5.0 (https://github.com/TransDecoder/TransDecoder)  
   MAGUS (https://github.com/vlasmirnov/MAGUS)  
   MAFFT v7.407 (https://mafft.cbrc.jp/alignment/software/)  
   ClipKIT (https://github.com/JLSteenwyk/ClipKIT)  
   trimAl v1.4.1 (http://trimal.cgenomics.org/)  
   BMGE v1.12 (https://anaconda.org/bioconda/bmge/files)  
   IQ-TREE v2.1.3 (https://github.com/iqtree/iqtree2)  
   ASTRAL-III v5.6.1 (https://github.com/smirarab/ASTRAL)  
   paml4.9j (http://abacus.gene.ucl.ac.uk/software/)  
   GNU Parallel 2018 (http://www.gnu.org/software/parallel/)

## User manual

The requirements for each script have been described in the beginning of the bash script text. Values of variables/parameters, such as tool path, number of threads etc., can be modified at the beginning of analyses following the script guidance.

Details of all commands and scripts usage are shown below. Operating system is popular Linux ones with ‘BASH’ shell, such as Centos, UBUNTU etc.

 ● BUSCO_extraction.sh:

1. Type 'bash BUSCO_extraction.sh BUSCO_folder', e.g., bash BUSCO_extraction.sh BUSCOs/. All the BUSCO results (run_* folders) are deposited in the same folder, e.g., BUSCOs/.
2. Tools TransDecoder, parallel are used in this script and will be automatically checked prior to formal analyses.
3. Modifying the head name of the fasta files for each locus, and placed all them in 0-raw_busco/.
4. Merging sequences of the same locus into the fasta files, and placed all these fasta files in 1-raw_loci/.
5. Filtering loci having too few taxa (less than three), and palced the rest of the fasta files in 2-loci_filter/.

 ● align_MAFFT.sh:

1. Type 'bash align_MAFFT.sh', e.g., bash align_MAFFT.sh.
2. Tools MAFFT and MAGUS are used in this script and will be automatically checked prior to formal analyses.
3. Input the option for MAFFT-based strategy: 1. mafft-auto; 2. linsi; 3. einsi; 4. ginsi; 5. MAGUS. Enter the number which strategy will be choosen.
4. Input the number of threads/cores (e.g., 8).
5. Input the name of input directory, e.g., 2-loci_filter/faa/.
6. Input the name of output directory, e.g., 3-faa_align. All the fasta files after aligning will be palced in this folder.

 ● trmming_alignments.sh:
 
1. Type 'bash trmming_alignments.sh', e.g., bash trmming_alignments.sh.
2. Tools trimAl, BMGE and ClipKIT are used in this script and will be automatically checked prior to formal analyses. 
3. Input the number of threads/cores (e.g., 8).
4. Input the option for input alignments: 1. amino acid; 2. nucleotide. Enter the number which alignments will be choosen.
5. Input the option for trimming tool: 1. trimAl; 2. BMGE; 3. ClipKIT. Enter the number which trimming tools will be choosen.
6. After selecting the appropriate tool, enter the corresponding file and the absolute path as prompted.
7. Input the name of output directory, e.g., 4-trim. All the fasta files after trimming will be palced in this folder.

 ● loci_filtering_alignment-based.sh:
 
1. Type 'bash loci_filtering_alignment-based.sh', e.g., bash loci_filtering_alignment-based.sh.
2. Tools PhyKIT and IQ-TREE are used in this script and will be automatically checked prior to formal analyses.
3. Input the number of threads/cores (e.g., 8).
4. Input the name of input directory containing all alignments, e.g., 4-trim/clipkit-kpic/.
5. Input the name of output directory, or an existing directory, e.g., length (the name of output directory can write according to the alignment-based strategy).
6. Input the option for alignment-based strategy for loci filtering: 1. alignment length; 2. number of parsimony-informative sites; 3. percentage of parsimony-informative sites in the alignment; 4. GC content; 5. RCV (Relative composition variability); 6. evolutionary rate (average pairwise identity); 7. likelihood mapping; 8. symmetry tests against SRH hypotheses. Enter the number which alignment-based strategy for loci filtering will be choosen.
7. Input the threshod. 





## Contact

Please send emails to Dr. Feng Zhang (xtmtd.zf@gmail.com).
