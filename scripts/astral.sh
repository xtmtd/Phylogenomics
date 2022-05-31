#!/bin/bash
#2021.08.05 by ZF

#prepare a list containing all loci names or gene tree names and a folder containing all gene tree files 
#ASTRAL is required
#generate the species tree usig ASTRAL
#Type "bash astral.sh"


#check ASTRAL
until [ -s $DIR_ASTRAL/astral*.jar ]
    do
      read -p "ASTRAL executable file is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/ASTRAL-5.7.1):      " DIR_ASTRAL_TEMP
      DIR_ASTRAL=$(realpath $(echo $DIR_ASTRAL_TEMP | sed "s/'//g"))
    done
echo "ASTRAL ...... OK"


#input the location of a list containing all loci names or gene tree names 
read -p "Please input a list containing all loci names or gene tree names, one name one line, e.g. /PATH/loci.ABS:      " LOCI_LIST_TEMP
LOCI_LIST=$(realpath $(echo $LOCI_LIST_TEMP | sed "s/'//g"))


#input the name of input directory
read -p "Please input the name of input directory containing all gene trees, e.g. /PATH/gene_trees:      " DIR_INPUT_TEMP
DIR_INPUT=$(realpath $(echo $DIR_INPUT_TEMP | sed "s/'//g"))


#input the name of output directory
read -p "Please input the name of output directory, or an existing directory:      " DIR_OUTPUT_TEMP
DIR_OUTPUT_TEMP1=$(echo $DIR_OUTPUT_TEMP | sed "s/'//g")
test -d $DIR_OUTPUT_TEMP1 || mkdir -p $DIR_OUTPUT_TEMP1
cd $DIR_OUTPUT_TEMP1 && DIR_OUTPUT=$(echo $PWD) && cd $DIR_OUTPUT


#merge all the gene trees into a single file
for loci in $(cat $LOCI_LIST); do cat $DIR_INPUT/"$loci"* >> all.gene.tre; done


#infer the ASTRAL tree
echo
echo "Calculating Astral tree ......"
java -jar $DIR_ASTRAL/astral*jar -i all.gene.tre -o species_tree.tre 2> log.txt


echo -e '\n'
echo "The file 'species_tree.tre' in the output directory is the MSC tree inferred from ASTRAL."


