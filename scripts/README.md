
# Phologenomics

## Introduction



All nine bash scripts are tested in Linux systems. Their functions are described below.

1. BUSCO_extraction.sh: this script extracts the amino acid and nucleotides sequences for all the BUSCO results (run_* folders) by using TransDecoder, i.e., modifying the head name of the fasta files for each locus and merging sequences of the same locus into the fasta files, then filtering loci having too few taxa (less than three).
2. align_MAFFT.sh: this script performs multiple sequence alignment by using MAFFT or MAGUS.
3. trmming_alignments.sh: this script performs trimming of multiple sequence alignments by using trimAl, BMGE or ClipKIT, i.e., identifying and removing highly divergent sites, or retaining parsimony-informative/constant sites from multiple sequence alignments.
4. loci_filtering_alignment-based.sh: this script filters loci by detecting alignment length, number of parsimony-informative sites, percentage of parsimony-informative sites in the alignment, GC content, compositional heterogeneity (RCV, Relative composition variability), evolutionary rate (average pairwise identity) and likelihood mapping criteria by using PhyKIT or IQ-TREE.
5. gene_trees.sh: this script generates individual gene trees for a set of alignments by using IQ-TREE.
6. loci_filtering_tree-based.sh: this script filter loci using gene tree-based methods, i.e., average bootstrap support (ABS), Degree of violation of the molecular clock (DVMC), treeness, signal-to-noise ratio (treeness over rcv), and spurious homologs (possible paralogs, incorrectly assembled sequences) by using PhyKIT or IQ-TREE.
7. matrix_generation.sh: this script generates the supermaxtrix, partition and occupancy for alignments by using PhyKIT.
8. astral.sh: this script generates the species tree by using ASTRAL.
9. mcmctree_AA.sh: this script estimates divergence time by using MCMCTree for large amino acid dataset.

## Requirements

Some bioinformatic tools are neccessary for above scripts. Most of them are recommended to be added into the environmental paths. Softwares, versions and source are listed below. All the softwares must be installed in Linux systems.

BBTools v38.32 (https://sourceforge.net/projects/bbmap/) 

Minia v3.2.4 (https://github.com/GATB/minia) 

Redundans v0.14c (https://github.com/lpryszcz/redundans) 

Minimap2 v2.12 (https://github.com/lh3/minimap2) 

Samtools v1.9 (http://www.htslib.org/) 

BESST v2.2.8 (https://github.com/ksahlin/BESST) 

GapCloser v1.12 (http://soap.genomics.org.cn/) 

BUSCO v3.0.2 (http://busco.ezlab.org/) 

FASconCAT-G v1.04 (https://github.com/PatrickKueck/FASconCAT-G) 

PHYLUCE v1.6.6 (http://phyluce.readthedocs.io/en/latest/index.html) 

faToTwoBit (http://hgdownload.soe.ucsc.edu/admin/exe/) 

ART-20160605 (https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) 

Stampy v1.0.32 (http://www.well.ox.ac.uk/project-stampy) 

PhyKIT (https://github.com/JLSteenwyk/PhyKIT)

R v4.1.2 (https://www.r-project.org/) 

TransDecoder v5.5.0 (https://github.com/TransDecoder/TransDecoder)

TreeShrink v1.3.8 (https://github.com/uym2/TreeShrink/tree/v1.3.8)

MAGUS (https://github.com/vlasmirnov/MAGUS) 

MAFFT v7.407 (https://mafft.cbrc.jp/alignment/software/) 

ClipKIT (https://github.com/JLSteenwyk/ClipKIT) 

trimAl v1.4.1 (http://trimal.cgenomics.org/) 

IQ-TREE v2.1.3 (https://github.com/iqtree/iqtree2) 

ASTRAL-III (https://github.com/smirarab/ASTRAL) 

pbmpi-1.8c (https://github.com/bayesiancook/pbmpi)

paml4.9j (http://abacus.gene.ucl.ac.uk/software/) 

## User manual


