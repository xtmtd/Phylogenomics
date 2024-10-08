# matrix_generation_mtDNA.sh —— a script for constructing the mitogenetic matrices rapidly

## Introduction

Mitochondrial genome (mitogenome) is widely used in phylogenetic relationships at different taxonomic levels. This script for constructing the mitogenetic matrices rapidly, integrating a series of computationally efficient bioinformatic tools, and were used with a universal ‘BASH’ shell or visual interface by Windows-like ‘drag and drop’ operations in LINUX systems.

This procedure was completed in a few min on desktop PCs (about 20 species). It provides convenience for accurate and rapid construction of matrices, before reconstructing the phylogentic trees.


## Requirements

Some bioinformatic tools are neccessary for above scripts. Most of them are recommended to be added into the environmental paths. Softwares, versions and source ate listed below.

FASconCAT-G v1.05 (https://github.com/PatrickKueck/FASconCAT-G)

MAFFT v7.520 (https://mafft.cbrc.jp/alignment/software)

trimAl v1.4.1 (https://github.com/inab/trimal)

Seqkit v2.4.0 (https://github.com/shenwei356/seqkit)

csvtk v0.25.0 (https://github.com/shenwei356/csvtk)

parallel (http://www.gnu.org/software/parallel)

TransDecoder v5.7.0 (https://github.com/TransDecoder/TransDecoder)

PhyKIT v1.11.15 (https://github.com/JLSteenwyk/PhyKIT)

IQ-TREE v2.2.2.7 (https://github.com/iqtree/iqtree2)


## User manual

Copy all the MitoZ annotation results, i.e. folders *.result, into the same folder (e.g. Mitoz_results).

Type 'bash matrix_generation_mtDNA.sh Mitoz_results'


## Contact

Please send emails to Dr. Feng Zhang (xtmtd.zf@gmail.com).
