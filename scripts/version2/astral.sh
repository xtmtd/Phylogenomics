#!/bin/bash
#2021.08.05 by ZF
#2023.06.06 by DSY
#2023.08.25 by DSY

#prepare a list containing all loci names or gene tree names and a folder containing all gene tree files 
#ASTER is required
#Type "bash astral.sh"


#input the name of output directory
read -p "Please input the name of output directory, or an existing directory:      " DIR_OUTPUT_TEMP
DIR_OUTPUT_TEMP1=$(echo ${DIR_OUTPUT_TEMP} | sed "s/'//g")
test -d ${DIR_OUTPUT_TEMP1} || mkdir -p ${DIR_OUTPUT_TEMP1}
cd ${DIR_OUTPUT_TEMP1} && DIR_OUTPUT=$(echo ${PWD}) && cd ${DIR_OUTPUT}


#Check the astral methods
read -p "Please input the option the strategy for astral: 1. ASTRAL; 2. weighted ASTRAL (wASTRAL); 3. CASTER-site; 4. CASTER-pair     " ASTER_METHOD
  until [ ${ASTER_METHOD} -gt 0 -a ${ASTER_METHOD} -lt 5 ]
    do
      read -p "Please input the option the strategy for astral: 1. ASTRAL; 2. weighted ASTRAL (wASTRAL); 3. CASTER-site; 4. CASTER-pair     " ASTER_METHOD
    done

#Check the threads can be used
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
  until [ ${THREADS} -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done


if [[ ${ASTER_METHOD} -eq 1 || ${ASTER_METHOD} -eq 2 ]]; then
  #input the name of input directory for gene trees
  echo
  read -p "Please input the name of input directory containing all gene trees, e.g. gene_trees:      " DIR_INPUT_TREE_TEMP
  DIR_INPUT_TREE=$(realpath $(echo ${DIR_INPUT_TREE_TEMP} | sed "s/'//g"))
  until [ -d ${DIR_INPUT_TREE} ]
     do
       read -p "Please input the name of input directory containing all gene trees, e.g. gene_trees:      " DIR_INPUT_TREE_TEMP
       DIR_INPUT_TREE=$(realpath $(echo ${DIR_INPUT_TREE_TEMP} | sed "s/'//g"))
     done
  echo "input directory containing gene trees ...... OK"


  #input the location of a list containing all loci names or gene tree names 
  read -p "Please input a list containing all loci names, one name one line, e.g. /PATH/loci.ABS:      " LOCI_LIST_TEMP
  LOCI_LIST=$(realpath $(echo ${LOCI_LIST_TEMP} | sed "s/'//g"))
  until [ -s ${LOCI_LIST} ]
     do
       read -p "Please input a list containing all loci names, one name one line, e.g. /PATH/loci.ABS:      " LOCI_LIST_TEMP
       LOCI_LIST=$(realpath $(echo ${LOCI_LIST_TEMP} | sed "s/'//g"))
     done
  echo "a list containing all loci names ...... OK"

  #merge all the gene trees into a single file
  for loci in $(cat ${LOCI_LIST}); do cat ${DIR_INPUT_TREE}/${loci}* >> all.gene.tre; done
  echo
  echo "Calculating ASTRAL/wASTRAL tree ......"
  echo

  #ASTRAL
  if [ ${ASTER_METHOD} -eq 1 ]; then
    #check ASTRAL
    if [ $(which astral) ]
      then
        echo "astral ...... OK"
        EXE_ASTER=$(which astral)
        DIR_ASTER=${EXE_ASTER%/*}
      else
        until [ -x ${DIR_ASTER}/astral ]
          do
            read -p "astral is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/ASTER/bin):      " DIR_ASTER_TEMP
          DIR_ASTER=$(echo ${DIR_ASTERAL_TEMP} | sed "s/'//g")
          done
        echo "astral ...... OK"
    fi

    #infer the ASTRAL tree
    astral -t ${THREADS} -r 16 -s 16 -i all.gene.tre -o species_tree.tre 2> log.txt

  #WASTRAL
  elif [ ${ASTER_METHOD} -eq 2 ]; then
    #check wASTRAL
    if [ $(which astral-hybrid) ]
      then
        echo "astral-hybrid ...... OK"
        EXE_WASTER=$(which astral-hybrid)
        DIR_WASTER=${EXE_WASTER%/*}
      else
        until [ -x ${DIR_WASTER}/astral-hybrid ]
          do
            read -p "astral-hybrid is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/ASTER-1.15/bin):      " DIR_WASTER_TEMP
            DIR_WASTER=$(echo ${DIR_WASTER_TEMP} | sed "s/'//g")
          done
        echo "astral-hybrid ...... OK"
    fi

    ${DIR_WASTER}/astral-hybrid -t ${THREADS} -r 16 -s 16 -x 100 -n 33 -i all.gene.tre -o species_tree.tre 2> log.txt

  fi


elif [[ ${ASTER_METHOD} -eq 3 || ${ASTER_METHOD} -eq 4 ]]; then
  #input the fasta file contaning all sequences
  read -p "Please input the fasta file contaning all sequences, e.g. /home/zf/Desktop/materials/test/matrix.fa:      " INPUT_FASTA_TEMP
  INPUT_FASTA=$(realpath $(echo ${INPUT_FASTA_TEMP} | sed "s/'//g"))
  until [ -s ${INPUT_FASTA} ]
     do
       read -p "Please input the fasta file contaning all sequences, e.g. /home/zf/Desktop/materials/test/matrix.fa:      " INPUT_FASTA_TEMP
       INPUT_FASTA=$(realpath $(echo ${INPUT_FASTA_TEMP} | sed "s/'//g"))
     done
  echo "FASTA matrix ...... OK"

  echo
  echo "Calculating caster-site/caster-pair tree ......"
  echo

  #caster-site
  if [ ${ASTER_METHOD} -eq 3 ]; then
    #check caster-site
    if [ $(which caster-site) ]
      then
        echo "caster-site ...... OK"
        EXE_CASTERSITE=$(which caster-site)
        DIR_CASTERSITE=${EXE_CASTERSITE%/*}
      else
        until [ -x ${DIR_CASTERSITE}/caster-site ]
          do
            read -p "caster-site is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/ASTER-1.15/bin):      " DIR_CASTERSITE_TEMP
            DIR_CASTERSITE=$(echo ${DIR_CASTERSITE_TEMP} | sed "s/'//g")
          done
        echo "caster-site ...... OK"
    fi

    #infer the caster-site tree
    ${DIR_CASTERSITE}/caster-site -t ${THREADS} -r 16 -s 16 -i ${INPUT_FASTA} -o species_tree.tre 2> log.txt

  elif [ ${ASTER_METHOD} -eq 4 ]; then
  #check caster-pair
    if [ $(which caster-pair) ]
      then
        echo "caster-pair ...... OK"
        EXE_CASTERPAIR=$(which caster-pair)
        DIR_CASTERPAIR=${EXE_CASTERPAIR%/*}
      else
        until [ -x ${DIR_CASTERPAIR}/caster-pair ]
          do
            read -p "caster-pair is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/ASTER-1.15/bin):      " DIR_CASTERPAIR_TEMP
            DIR_CASTERPAIR=$(echo ${DIR_CASTERPAIR_TEMP} | sed "s/'//g")
          done
        echo "caster-pair ...... OK"
    fi

    #infer the caster-pair tree
    ${DIR_CASTERPAIR}/caster-pair -t ${THREADS} -r 16 -s 16 -i ${INPUT_FASTA} -o species_tree.tre 2> log.txt

  fi

fi

echo
echo "The file 'species_tree.tre' in the output directory is the MSC tree inferred from ASTER."
