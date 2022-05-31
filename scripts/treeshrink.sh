#!/bin/bash
#2021.10.10 by DSY

#detection of outlier long branches in collections of phylogenetic trees
#type "bash treeshrink.sh"
#treeshrink, python and R is required



DIR_CURR=$(echo $PWD)

#input the name of input directory containing all gene trees
read -p "Please input the name of input directory containing all gene trees, e.g. gene_trees:      " DIR_INPUT_TREE_TEMP
DIR_INPUT_TREE_TEMP1=$(echo $DIR_INPUT_TREE_TEMP | sed "s/'//g")
cd $DIR_INPUT_TREE_TEMP1 && DIR_INPUT_TREE=$(echo $PWD) && cd $DIR_CURR


#input the name of input directory containing all loci alignments
read -p "Please input the name of input directory containing all loci alignments, e.g. matrix75/alignment:      " DIR_INPUT_ALIGN_TEMP
DIR_INPUT_ALIGN_TEMP1=$(echo $DIR_INPUT_ALIGN_TEMP | sed "s/'//g")
cd $DIR_INPUT_ALIGN_TEMP1 && DIR_INPUT_ALIGN=$(echo $PWD) && cd $DIR_CURR


#input the number of α threshold (α = 0.05 is default)
read -p "Please input the number of α threshold (e.g. 0.05):      " THRESHOLD
  until [ $THRESHOLD ]
    do
      read -p "Please input the number of α threshold (e.g. 0.05):      " THRESHOLD
    done

#copy the alignments and gene trees, and extract the loci
mkdir -p treeshrink-fas treeshrink/treeshrink_fas treeshrink/treeshrink_treefile treeshrink/shrunk$THRESHOLD treeshrink-treefile
cp $DIR_INPUT_ALIGN/* treeshrink-fas
cp $DIR_INPUT_TREE/* treeshrink-treefile
ls treeshrink-fas | sed 's/.fas//g' > treeshrink/loci.list

LOCI_NAME=$(cat treeshrink/loci.list)
for loci in $LOCI_NAME
do 
  mkdir treeshrink/shrunk$THRESHOLD/$loci/
  cp treeshrink-fas/$loci.fas treeshrink/shrunk$THRESHOLD/$loci/
  cp treeshrink-treefile/$loci.fas.treefile treeshrink/shrunk$THRESHOLD/$loci/
done

rm -rf treeshrink-treefile/ treeshrink-fas/

#run treeshrink
for loci in $LOCI_NAME
  do 
    run_treeshrink.py -i treeshrink/shrunk$THRESHOLD/ -t $loci.fas.treefile -q "$THRESHOLD" -a $loci.fas -O $loci.treeshrink
    mv treeshrink/shrunk$THRESHOLD/$loci.treeshrink_summary.txt treeshrink/shrunk$THRESHOLD/$loci/
    cp treeshrink/shrunk$THRESHOLD/$loci/$loci.treeshrink.fas treeshrink/treeshrink_fas
    cp treeshrink/shrunk$THRESHOLD/$loci/$loci.treeshrink.treefile treeshrink/treeshrink_treefile
  done