# Scripts of phylogenomics

## Introduction

We present a set of custom scripts for aligning, trimming, filtering and creating concatenation matrix before reconstructing the phylogenetic trees, employing a series of computationally efficient bioinformatic tools. All the scripts must be used in LINUX systems. The analysis process is shown in the following figure. These scripts have major advantages in data analysis automatically, and reducing computing consumption. Most scripts are already parallelized. It provides an important process for phylogenomics.

![image](https://user-images.githubusercontent.com/45136134/199656510-200f05e3-2227-428a-97ac-69f0b3d96fb6.png)

All nine bash scripts are tested in Linux systems, such as CentOS 7, CentOS 8, ubuntu18.04 LTS, ubuntu20.04 LTS, ubuntu22.04 LTS etc. Their functions are described below.

1. **BUSCO_extraction.sh:** this script extracts the amino acid and nucleotides sequences for all the BUSCO results (run_* folders) by using TransDecoder, i.e., modifying the head name of the fasta files for each locus and merging sequences of the same locus into the fasta files, then filtering loci having too few taxa (less than three).
2. **align_MAFFT.sh:** this script performs multiple sequence alignment by using MAFFT or MAGUS.
3. **trimming_alignments.sh:** this script performs trimming of multiple sequence alignments by using trimAl, BMGE or ClipKIT, i.e., identifying and removing highly divergent sites, or retaining parsimony-informative/constant sites from multiple sequence alignments.
4. **loci_filtering_alignment-based.sh:** this script filters loci by detecting alignment length, number of parsimony-informative sites, percentage of parsimony-informative sites in the alignment, GC content, compositional heterogeneity (RCV, Relative composition variability), evolutionary rate (average pairwise identity) and likelihood mapping criteria by using PhyKIT or IQ-TREE.
5. **gene_trees.sh:** this script generates individual gene trees for a set of alignments by using IQ-TREE.
6. **loci_filtering_tree-based.sh:** this script filter loci using gene tree-based methods, i.e., average bootstrap support (ABS), Degree of violation of the molecular clock (DVMC), treeness, signal-to-noise ratio (treeness over rcv), spurious homologs (possible paralogs, incorrectly assembled sequences), and treeshrink by using PhyKIT, IQ-TREE or TreeShrink.
7. **matrix_generation.sh:** this script generates the supermaxtrix, partition and occupancy for alignments by using PhyKIT.
8. **astral.sh:** this script generates the species tree by using ASTRAL.
9. **mcmctree_AA.sh:** this script estimates divergence time by using MCMCTree for "large scale" amino acid dataset.

## Requirements

Some bioinformatic tools are neccessary for above scripts. Most of them are recommended to be added into the environmental paths. Softwares, versions and source are listed below. All the softwares must be installed in LINUX systems.

   FASconCAT-G v1.04 (https://github.com/PatrickKueck/FASconCAT-G)   
   PHYLUCE v1.6.6 (http://phyluce.readthedocs.io/en/latest/index.html)   
   PhyKIT v1.11.10 (https://github.com/JLSteenwyk/PhyKIT)   
   TransDecoder v5.5.0 (https://github.com/TransDecoder/TransDecoder)   
   MAGUS v0.1.0b0 (https://github.com/vlasmirnov/MAGUS)   
   MAFFT v7.407 (https://mafft.cbrc.jp/alignment/software/)   
   ClipKIT v1.1.5 (https://github.com/JLSteenwyk/ClipKIT)   
   trimAl v1.4.1 (http://trimal.cgenomics.org/)   
   BMGE v1.12 (https://anaconda.org/bioconda/bmge/files)   
   IQ-TREE v2.1.3 (https://github.com/iqtree/iqtree2)   
   ASTRAL-III v5.6.1 (https://github.com/smirarab/ASTRAL)   
   paml4.9j (http://abacus.gene.ucl.ac.uk/software/)   
   GNU Parallel 2018 (http://www.gnu.org/software/parallel/)   
   R v4.1.2 (https://www.r-project.org/)   
   TreeShrink v1.3.8b (https://github.com/uym2/TreeShrink)   
   
   
To faciliate the use of our pipeline and simplify the installation of various packages, we also prepared the virtual mirror (.vmdk), including all the bioinformatic tools and custom scripts ready in the CentOS 7 system. It can be directly opened by VMware or VirtualBox, which are often installed on Windows systems. Please download it here: https://dx.doi.org/10.6084/m9.figshare.21283026. More than 100 Gb (at least 50) disk space is recommented for its use. Below is related information on CentOS system:

operating system: CentOS 7

passwords of accounts, i.e. root and zf: 1

Most packages are installed under /usr/local/bin, /home/zf/install/ etc. Environmental variables and paths can be checked by "vi ~/.bashrc" and "printenv". Most executables can be directly performed without typing the installation paths. UCE-related analyses should be executed by initiating conda environments (/home/zf/install/miniconda/envs/), such as phyluce, just type 'source activate phyluce'.

The video introduction and “walk-through” for each script are provided at https://space.bilibili.com/319699648/channel/seriesdetail?sid=2682055.

## User manual

The requirements for each script have been described in the beginning of the bash script text. Values of variables/parameters, such as tool path, number of threads etc., can be modified at the beginning of analyses following the script guidance.

Details of all scripts usage are shown below. Operating system is popular Linux ones with ‘BASH’ shell, such as Centos, UBUNTU etc.

** *NOTE*: Do not change the file name (loci name) arbitrarily! **

** *SUGGSET*: After each step of analysing, please check the number of loci! **

 ● **BUSCO_extraction.sh:**

1. Type 'bash BUSCO_extraction.sh'.
2. Tools TransDecoder, parallel are used in this script and will be automatically checked prior to formal analyses.
3. Input the number of threads/cores (e.g., 8).
4. Input the name of input directory containing all alignments, e.g., /PATH/BUSCOs/ (All the BUSCO results (run_* folders) are deposited in the same folder, e.g., BUSCOs/).
5. Modify the head name of the fasta files for each locus, and placed all them in 0-raw_busco/.
6. Merge sequences of the same locus into the fasta files, and placed all these fasta files in 1-raw_loci/.
7. Filter loci having too few taxa (less than three), and placed the rest of the fasta files in 2-loci_filter/.

 ● **align_MAFFT.sh:**

1. Type 'bash align_MAFFT.sh'.
2. Tools MAFFT and MAGUS are used in this script and will be automatically checked prior to formal analyses.
3. Input the option for MAFFT-based strategy: 1. mafft-auto; 2. linsi; 3. einsi; 4. ginsi; 5. MAGUS. Enter the number which strategy be choosen.
4. Input the number of threads/cores (e.g., 8).
5. Please input the option for input unalignments: 1. amino acid; 2. nucleotide. Enter the number which alignments be choosen.
6. Input the name of input directory containing all alignments, e.g., 2-loci_filter/faa/.
7. input the name of output directory, or an existing directory, e.g., 3-faa_align. All the fasta files after aligning will be placed in this folder.
8. Do you want to align the the corresponding nucleotide sequences? 1. Yes; 2. No. 
9. Please input the folder names (with its path) containing unaligned nucleotide sequences, for example '2-loci_filter/fna': 2-loci_filter/fna/. (This step is performed only when 1 [Yes] is selected in Step 8)

 ● **trimming_alignments.sh:**
 
1. Type 'bash trimming_alignments.sh'.
2. Tools parallel, trimAl, BMGE, ClipKIT and PhyKIT are used in this script and will be automatically checked prior to formal analyses. 
3. Input the number of threads/cores (e.g., 8).
4. Input the option for input alignments: 1. amino acid; 2. nucleotide. Enter the number which alignments be choosen.
5. Input the option for trimming tool: 1. trimAl; 2. BMGE; 3. ClipKIT. Enter the number which trimming tools be choosen.
6. After selecting the appropriate tool, enter the corresponding file and the absolute path as prompted.
7. Input the name of output directory, e.g., 4-trim/. All the fasta files after trimming will be placed in this folder.
8. Do you want to trim the the corresponding nucleotide sequences? 1. Yes; 2. No.
9. Please input the folder names (with its path) containing nucleotide alignments, for example '3-faa_align/fna', e.g., 3-faa_align/fna/. (This step is performed only when 1 [Yes] is selected in Step 8)

 ● **loci_filtering_alignment-based.sh:**
 
1. Type 'bash loci_filtering_alignment-based.sh'.
2. Tools parallel, PhyKIT and IQ-TREE are used in this script and will be automatically checked prior to formal analyses.
3. Input the number of threads/cores (e.g., 8).
4. Input the name of input directory containing all alignments, e.g., 4-trim/clipkit-kpic/.
5. Input the name of output directory, or an existing directory, e.g., length (the name of output directory can write according to the alignment-based strategy). All the fasta files after filtering based alignments, and the list of loci will be placed in this folder.
6. Input the option for alignment-based strategy for loci filtering: 1. alignment length; 2. number of parsimony-informative sites; 3. percentage of parsimony-informative sites in the alignment; 4. GC content; 5. RCV (Relative composition variability); 6. evolutionary rate (average pairwise identity); 7. likelihood mapping; 8. symmetry tests against SRH hypotheses. Enter the number which alignment-based strategy for loci filtering be choosen.
7. Input the threshold, i.e., the length or number of parsimony-informative sites threshold (the minimum value), the RCV threshold (the maximum value), the GC percentage threshold (the maximum value), the evolutionary rate threshold (the maximum value), the FcLM quartet threshold (the minimum value), the symtest p-value threshold (the maximum value).

 ● **gene_trees.sh:**
 
1. Type 'bash gene_trees.sh'.
2. Tools IQ-TREE and parallel are used in this script and will be automatically checked prior to formal analyses.
3. Input the name of input directory containing all alignments, e.g., 4-trim/clipkit-kpic/.
4. Input the name of output directory, or an existing directory, e.g., gene_trees. All gene trees have been deposited in this folder.
5. Input the option for input alignments: 1. amino acid; 2. nucleotide. Enter the number which alignments be choosen.
6. Input the option for protein substitution model: 1. automatically ModelFinder searching; 2. restricted to LG model to reduce computational burden (-mset LG); 3. mixture model EX_EHO (-mset EX_EHO, may be very time-consuming but more accurate); 4. restricted to LG/EX_EHO model to reduce computational burden (-mset LG,EX_EHO); 5. insect universal model Q.insect (may be better than LG model for insects); 6. mixture model LG+C20+FO+R (-mset LG+C20+FO+R, very time-comsuming); 7. mitochondrial (--msub mitochondrial). Or input the option for DNA substitution model: 1. automatically ModelFinder searching; 2. restricted to HKY/GTR model to reduce computational burden (-mset HKY,GTR). Enter the number which protein or DNA substitution model be choosen.
7. Input the number of threads/cores used for each IQ-TREE analysis (e.g., 1).
8. Input the number of IQ-TREE jobs/tasks (e.g., 8).

 ● **loci_filtering_tree-based.sh:**
 
1. Type 'bash loci_filtering_tree-based.sh'.
2. Tools parallel, R, TreeShrink, PhyKIT and Seqkit are used in this script and will be automatically checked prior to formal analyses.
3. Input the number of threads/cores (e.g., 8).
4. Input the name of input directory containing all gene trees, e.g., gene_trees/.
5. Input the name of input directory containing all loci alignments, e.g. 4-trim/clipkit-kpic/.
6. Input the name of output directory, or an existing directory, e.g., ABS (the name of output directory can write according to the tree-based strategy). All the fasta files after filtering based gene trees, and the list of loci will be placed in this folder.
7. Input the option for tree-based strategy for loci filtering: 1. average bootstrap support (ABS); 2. degree of violation of the molecular clock (DVMC); 3. treeness; 4. signal-to-noise ratio (treeness over RCV); 5. spurious homologs identification; 6. treeshrink. Enter the number which tree-based strategy for loci filtering be choosen.
8. Input the threshold, i.e., the ABS threshold (the minimum value), the DVMC threshold (the maximum value), the treeness or treeness/RCV threshold (the minimum value), the treeshrink α threshold.

 ● **matrix_generation.sh:**
 
1. Type 'bash matrix_generation.sh'.
2. Tools Seqkit and PhyKIT are used in this script and will be automatically checked prior to formal analyses.
3. Input the name of input directory containing all alignments, e.g., 4-trim/clipkit-kpic/.
4. Input the name of output directory, or an existing directory, e.g., matrix90.
5. Input the name of a prefix for the generated matrix-related files, e.g., DATASET1.
6. Input the minimum percentage value for taxa occupancy, usually ranging from 50% to 100% (e.g., 50, 75, 90): 90.
7. Move the outgroup species to the first one in the alignment matrix? 1. Yes; 2. No (Many phylogenetic tools, such as IQ-TREE, will view the first taxon of the alignment as the 'outgroup' in the final tree file). If 'Yes' (i.e., input '1'), please input the species name of outgroup taxon, such as 'Zootermopsis_nevadensis'.
8. Individual loci alignments, concatenated matrix and partition file are deposited in the OUTPUT/matrix$OCCUPANCY.

 ● **astral.sh:**
 
1. Type 'bash astral.sh'.
2. Tool ASTRAL is used in this script and will be automatically checked prior to formal analyses.
3. Input a list containing all loci names or gene tree names, one name one line, e.g., /PATH/loci.ABS.
4. Input the name of input directory containing all gene trees, e.g., /PATH/gene_trees/.
5. Input the name of output directory, or an existing directory, e.g., ASTRAL.
6. The file 'species_tree.tre' in the output directory is the MSC tree inferred from ASTRAL.

 ● **mcmctree_AA.sh:**
 
1. Type 'bash mcmctree_AA.sh'.
2. Tools PhyKIT, PAML, trimAl, parallel and csvtk are used in this script and will be automatically checked prior to formal analyses.
3. Input the number of threads/cores (e.g., 8).
4. Default burnin and sample frequency (sampfreq) are set as 50000 and 5 in this script. Please input the number of samples/generations kept for MCMCTree, 10000~20000 could be enough for most cases), e.g., 20000.
5. Input the name of input directory containing all alignments, e.g., /PATH/alignments.
6. Input the partition file 'XXXbest_model.nex' (with its path) generated by IQ-TREE, e.g., /PATH/test.partition.best_model.nex.


## Contact

Please send emails to Dr. Shiyu Du (zjjhdsy@126.com) or Prof. Feng Zhang (xtmtd.zf@gmail.com).
