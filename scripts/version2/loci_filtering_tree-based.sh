#!/bin/bash
#2021.08.04 wrote by ZF
#2023.07.13 revised by DSY
#2023.08.25 by DSY
#2023.10.22 by ZF

#filter loci using gene tree-based methods, i.e. average bootstrap support (ABS), Degree of violation of the molecular clock (DVMC), treeness, signal-to-noise ratio (treeness over rcv), spurious homologs (possible paralogs, incorrectly assembled sequences), treeshrink, and inconsistent genes
#Type 'bash loci_filtering_tree-based.sh'
#parallel, cscvtk, TreeShrink, and PhyKIT may be required


#check Phykit
echo
if [ $(which phykit) ]
    then
      echo "PhyKIT ...... OK"
      EXE_PHYKIT=$(which phykit)
      DIR_PHYKIT=${EXE_PHYKIT%/*}
    else
      until [ -x ${DIR_PHYKIT}/phykit ]
        do
          read -p "PhyKIT is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_PHYKIT_TEMP
          DIR_PHYKIT=$(realpath $(echo ${DIR_PHYKIT_TEMP} | sed "s/'//g"))
        done
      echo "PhyKIT ...... OK"
fi


#check parallel
echo
if [ $(which parallel) ]
    then
      echo "parallel ...... OK"
      EXE_PARALLEL=$(which parallel)
      DIR_PARALLEL=${EXE_PARALLEL%/*}
    else
      until [ -x ${DIR_PARALLEL}/parallel ]
        do
          read -p "parallel is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/TransDecoder-v5.5.0):      " DIR_PARALLEL_TEMP
          DIR_PARALLEL=$(realpath $(echo ${DIR_PARALLEL_TEMP} | sed "s/'//g"))
        done
      echo "PARALLEL ...... OK"
fi


#check csvtk
echo
if [ $(which csvtk) ]
    then
      echo "csvtk ...... OK"
      EXE_CSVTK=$(which csvtk)
      DIR_CSVTK=${EXE_CSVTK%/*}
    else
      until [ -x ${DIR_CSVTK}/csvtk ]
  do
    read -p "CSVTK is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_CSVTK_TEMP
    DIR_CSVTK=$(realpath $(echo ${DIR_CSVTK_TEMP} | sed "s/'//g"))
  done
      echo "csvtk ...... OK"
fi


#Check the threads can be used
echo
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
  until [ ${THREADS} -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done

#input the name of input directory for gene trees
echo
read -p "Please input the name of input directory containing all gene trees, e.g. gene_trees:      " DIR_INPUT_TREE_TEMP
DIR_INPUT_TREE=$(realpath $(echo ${DIR_INPUT_TREE_TEMP} | sed "s/'//g"))
  until [ -d ${DIR_INPUT_TREE} ]
     do
       read -p "Please input the name of input directory containing all gene trees, e.g. gene_trees:      " DIR_INPUT_TREE_TEMP
       DIR_INPUT_TREE=$(realpath $(echo ${DIR_INPUT_TREE_TEMP} | sed "s/'//g"))
     done
  echo "input directory for gene trees ...... OK"


#input the name of input directory for all alignments
echo
read -p "Please input the name of input directory containing all loci alignments, e.g. 4-trim/clipkit-kpic:      " DIR_INPUT_ALIGN_TEMP
DIR_INPUT_ALIGN=$(realpath $(echo ${DIR_INPUT_ALIGN_TEMP} | sed "s/'//g"))
  until [ -d ${DIR_INPUT_ALIGN} ]
     do
       read -p "Please input the name of input directory containing all loci alignments, e.g. 4-trim/clipkit-kpic:      " DIR_INPUT_ALIGN_TEMP
       DIR_INPUT_ALIGN=$(realpath $(echo ${DIR_INPUT_ALIGN_TEMP} | sed "s/'//g"))
     done
  echo "input directory for all alignments ...... OK"



#input the name of output directory
read -p "Please input the name of output directory, or an existing directory:      " DIR_OUTPUT_TEMP
DIR_OUTPUT_TEMP1=$(echo ${DIR_OUTPUT_TEMP} | sed "s/'//g")
test -d ${DIR_OUTPUT_TEMP1} || mkdir -p ${DIR_OUTPUT_TEMP1}
cd ${DIR_OUTPUT_TEMP1} && DIR_OUTPUT=$(echo ${PWD}) && cd ${DIR_OUTPUT}


#Check the loci filtering method can be used
read -p "Please input the option for tree-based strategy for loci filtering: 1. average bipartition support (ABS); 2. degree of violation of the molecular clock (DVMC); 3. treeness; 4. signal-to-noise ratio (treeness over RCV); 5. evolutionary rate; 6. saturation; 7. spurious homologs/outliers identification; 8. inconsistent genes     " FILTER_METHOD
  until [[ $FILTER_METHOD -gt 0 && ${FILTER_METHOD} -lt 9 ]]
    do
      read -p "Please input the option for tree-based strategy for loci filtering: 1. average bipartition support (ABS); 2. degree of violation of the molecular clock (DVMC); 3. treeness; 4. signal-to-noise ratio (treeness over RCV);  5. evolutionary rate; 6. saturation; 7. spurious homologs/outliers identification; 8. inconsistent genes     " FILTER_METHOD
    done


#generate output folders and loci list
mkdir -p ${DIR_OUTPUT}/alignments_remaining ${DIR_OUTPUT}/gene_trees_remaining ${DIR_OUTPUT}/temp
cd ${DIR_INPUT_ALIGN} && ls * > ${DIR_OUTPUT}/loci.list #generate loci list
cd ${DIR_OUTPUT}


#average bipartition support (ABS)
if [ ${FILTER_METHOD} -eq 1 ]; then
  #calculate ABS values for gene trees
  echo
  echo "Calculating average bipartition support (ABS) ......"

  export DIR_INPUT_TREE DIR_PHYKIT
  temp_fun() {
    echo -e "$1""\t""$(${DIR_PHYKIT}/phykit bss ${DIR_INPUT_TREE}/"$1"* | awk 'NR==2 {print $2}')" > ./temp/$1
  }
  export -f temp_fun
  cat loci.list | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 temp_fun %
  cat temp/* | sort -k2nr > ABS.txt
  rm -rf temp

  #generate the frequency distribution plot and help to determine the threshold
  awk 'NF>1 && $2 != "" {print $2}' ABS.txt | 
sort -n | 
awk -v OFS="\t" 'BEGIN{print "Loci_number\taverage_bipartition_support"} {print NR, $0}' | 
csvtk -t plot line -x Loci_number -y "average_bipartition_support" --format pdf > ABS.pdf

  #input the threshold 
  echo
  echo "Read the file ABS.txt in the output folder to determine the minimum ABS threshold (usually 50~90). Loci of higher ABS values are thought to be of more phylogenetic signal. The distribution plot 'ABS.pdf' may also help to determine the cutoff values."
  echo
  read -p "Please input the threshold (the minimum value):      " THRESHOLD

  #generate the list of remaining loci
  awk '$2>='"${THRESHOLD}"' {print $1}' ABS.txt > loci.remaining


#degree of violation of the molecular clock (DVMC)
elif [ ${FILTER_METHOD} -eq 2 ]; then
  #input the file name contaning outgroup taxa names
  read -p "Please input the file name (with its complete path) contaning outgroup taxa names, one name per row, e.g. /home/zf/Desktop/materials/test/6-gene_trees/DVMC/root.txt:      " INPUT_ROOT_TEMP
  INPUT_ROOT=$(realpath $(echo ${INPUT_ROOT_TEMP} | sed "s/'//g"))

  #calculate DVMC values for gene trees
  echo
  echo "Calculating degree of violation of the molecular clock (DVMC) ......"

  export DIR_INPUT_TREE DIR_PHYKIT INPUT_ROOT
  temp_fun() {
    phykit tl ${DIR_INPUT_TREE}/"$1"* | grep -f ${INPUT_ROOT} > root."$1".txt    #get outgroup list of the gene tree
    test -s root."$1".txt && cp ${DIR_INPUT_TREE}/"$1"* "$1".pruned || (phykit prune ${DIR_INPUT_TREE}/"$1"* root."$1".txt -o "$1".pruned)    #prune tree
    echo -e "$1""\t""$(phykit dvmc "$1".pruned)" > ./temp/$1    #DVMC
  }
  export -f temp_fun
  cat loci.list | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 temp_fun %
  rm *.txt *.pruned
  cat temp/* | sort -k2n > DVMC.txt
  rm -rf temp

  #generate the frequency distribution plot and help to determine the threshold
  awk '{print $2}' DVMC.txt | sort -k1n | awk -v OFS="\t" '{print NR, $0}' | sed '1i\Loci_number\tDVMC' | csvtk -t plot line -x Loci_number -y "DVMC" --format pdf > DVMC.pdf

  #input the threshold 
  echo
  echo "Read the file DVMC.txt in the output folder to determine the maximum DVMC threshold. Lower DVMC values are thought to be desirable because they are indicative of a lower degree of violation in the molecular clock assumption. The distribution plot 'DVMC.pdf' may also help to determine the cutoff values."
  echo
  read -p "Please input the threshold (the maximum value):      " THRESHOLD

  #generate the list of remaining loci
  awk '$2<='"${THRESHOLD}"' {print $1}' DVMC.txt > loci.remaining


#treeness and treeness/rcv
elif [[ ${FILTER_METHOD} -eq 3 || ${FILTER_METHOD} -eq 4 ]]; then 
  #calculate basic values for gene trees
  echo
  echo "Calculating treeness and treeness/rcv ......"

  export DIR_INPUT_TREE DIR_INPUT_ALIGN DIR_PHYKIT
  temp_fun() {
    echo -e "$1""\t""$(${DIR_PHYKIT}/phykit tor -a ${DIR_INPUT_ALIGN}/$1 -t ${DIR_INPUT_TREE}/"$1"*)" > ./temp/$1
  }
  export -f temp_fun
  cat loci.list | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 temp_fun %
  cat temp/* | sort -k2nr > TOR.txt
  rm -rf temp

  #generate the frequency distribution plot and help to determine the threshold
  if [ ${FILTER_METHOD} -eq 3 ]; then
    awk '{print $3}' TOR.txt | sort -k1n | awk -v OFS="\t" '{print NR, $0}' | sed '1i\Loci_number\ttreeness' | csvtk -t plot line -x Loci_number -y "treeness" --format pdf > treeness.pdf
  elif [ ${FILTER_METHOD} -eq 4 ]; then
    awk '{print $2}' TOR.txt | sort -k1n | awk -v OFS="\t" '{print NR, $0}' | sed '1i\Loci_number\ttreeness over RCV' | csvtk -t plot line -x Loci_number -y "treeness over RCV" --format pdf > treeness_over_RCV.pdf
  else
    echo
  fi

  #input the threshold 
  echo
  echo "Read the file TOR.txt in the output folder to determine the minimum threshold. The values of 2nd~4th columns represent treeness/RCV, treeness, and RCV, respectively. Higher treeness/RCV or treeness values are thought to be desirable. The distribution plot 'XXX.pdf' may also help to determine the cutoff values."
  echo
  read -p "Please input the threshold (the minimum value):      " THRESHOLD

  #generate the list of remaining loci
  if [ ${FILTER_METHOD} -eq 3 ]; then
    awk '$3>='"${THRESHOLD}"' {print $1}' TOR.txt > loci.remaining
  elif [ ${FILTER_METHOD} -eq 4 ]; then
    awk '$2>='"${THRESHOLD}"' {print $1}' TOR.txt > loci.remaining
  else
    echo
  fi


#evolutionary rate
elif [ ${FILTER_METHOD} -eq 5 ]; then
  #calculate evolutionary rate values for gene trees
  echo
  echo "Calculating evolutionary rates (evo_rate) ......"

  export DIR_INPUT_TREE DIR_PHYKIT
  temp_fun() {
    echo -e "$1""\t""$(${DIR_PHYKIT}/phykit evo_rate ${DIR_INPUT_TREE}/"$1"*)" > ./temp/$1
  }
  export -f temp_fun
  cat loci.list | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 temp_fun %
  cat temp/* | sort -k2nr > evo_rate.txt
  rm -rf temp

  #generate the frequency distribution plot and help to determine the threshold
  awk '{print $2}' evo_rate.txt | sort -k1n | awk -v OFS="\t" '{print NR, $0}' | sed '1i\Loci_number\tevolutionary rate' | csvtk -t plot line -x Loci_number -y "evolutionary rate" --format pdf > evo_rate.pdf

  #input the threshold 
  echo
  echo "Read the file evo_rate.txt in the output folder to determine the rate threshold. Very high values indicate that those 'fast' genes may be more suitable for shallow phylogenies. The distribution plot 'evo_rate.pdf' may also help to determine the cutoff values."
  echo
  read -p "Please input the minimum evo_rate threshold:      " MIN_evo_rate
  echo
  read -p "Please input the maximum evo_rate threshold:      " MAX_evo_rate

  #generate the list of remaining loci
  awk '$2>='"${MIN_evo_rate}"'&&$2<='"${MAX_evo_rate}"' {print $1}' evo_rate.txt > loci.remaining


#saturation
elif [ ${FILTER_METHOD} -eq 6 ]; then
  #calculate ABS values for gene trees
  echo
  echo "Calculating saturation values ......"

  export DIR_INPUT_TREE DIR_INPUT_ALIGN DIR_PHYKIT
  temp_fun() {
    echo -e "$1""\t""$(${DIR_PHYKIT}/phykit saturation -a ${DIR_INPUT_ALIGN}/$1 -t ${DIR_INPUT_TREE}/"$1"*)" > ./temp/$1
  }
  export -f temp_fun
  cat loci.list | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 temp_fun %
  cat temp/* | sort -k2nr > saturation.txt
  rm -rf temp

  #generate the frequency distribution plot and help to determine the threshold
  awk '{print $2}' saturation.txt | sort -k1n | awk -v OFS="\t" '{print NR, $0}' | sed '1i\Loci_number\tsaturation' | csvtk -t plot line -x Loci_number -y "saturation" --format pdf > saturation.pdf

  #input the threshold 
  echo
  echo "Read the file saturation.txt in the output folder to determine the minimum the value of saturation. Loci of higher values of saturation are thought to be desirable. The distribution plot 'saturation.pdf' may also help to determine the cutoff values."
  echo
  read -p "Please input the threshold (the minimum value):      " THRESHOLD

  #generate the list of remaining loci
  awk '$2>='"${THRESHOLD}"' {print $1}' saturation.txt > loci.remaining


#detect spurious homologs/outliers
elif [ ${FILTER_METHOD} -eq 7 ]; then
  #check treeshrink
  if [ $(which run_treeshrink.py) ]
    then
      echo "run_treeshrink.py ...... OK"
      EXE_TREESHRINK=$(which run_treeshrink.py)
      DIR_TREESHRINK=${EXE_TREESHRINK%/*}
    else
      until [ -x ${DIR_TREESHRINK}/run_treeshrink.py ]
        do
          read -p "run_treeshrink.py is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/TreeShrink-1.3.8b):      " DIR_TREESHRINK_TEMP
          DIR_TREESHRINK=$(realpath $(echo ${DIR_TREESHRINK_TEMP} | sed "s/'//g"))
        done
      echo "treeshrink ...... OK"
  fi

  #Treeshrink for gene trees
  #input the number of α threshold (α = 0.05 is default)
  read -p "Please input the number of α threshold. Default is 0.05 (e.g. 0.05):      " THRESHOLD
    until [ ${THRESHOLD} ]
      do
        read -p "Please input the number of α threshold. Default is 0.05 (e.g. 0.05):      " THRESHOLD
      done

  echo
  echo "Calculating spurious homologs/outliers ......"
 
  #generate folders for each locus
  for id in $(cat loci.list)
    do
      mkdir -p temp/${id}
      cp ${DIR_INPUT_ALIGN}/${id} temp/${id}/input.fasta
      cp ${DIR_INPUT_TREE}/${id}* temp/${id}/input.tree
    done
  
  #treeshrink
  ${DIR_TREESHRINK}/run_treeshrink.py -i temp/ -t input.tree -q "${THRESHOLD}" -a input.fasta
  for id in $(cat loci.list); do test -s temp/${id}/output.txt && echo ${id} >> list.outlier; done #problematic loci

  #copy corresponding alignments and gene trees
  for id in $(cat loci.list)
    do
      cp temp/${id}/output.fasta alignments_remaining/${id}
      cp temp/${id}/output.tree gene_trees_remaining/${id}.treefile
    done

  #remaining alignments for those outliers
  mkdir alignments_remaining_outlier
  for id in $(cat list.outlier); do cp alignments_remaining/${id} alignments_remaining_outlier/; done #these filtered alignments may be used to re-construct phylogenetic trees, which are possibly more accurate those pruned by treeshrink

  echo
  echo "The file loci.outlier contains those alignments in which TreeShrink has detected spurious homologs/outliers."
  echo
  echo "All the alignments and gene trees are deposited in the folder 'alignments_remaining' and 'gene_trees_remaining', respectively. Spurious homologs/outliers for the problematic loci have been removed from alignments and gene trees."
  echo
  echo "The folder 'alignments_remaining_outlier' only includes those filtered alignments. They may be used to re-construct phylogenetic trees, which are possibly more accurate those pruned by TreeShrink."


#inconsistent genes
elif [ ${FILTER_METHOD} -eq 8 ]; then
  #check ASTRAL
  until [ -s ${DIR_ASTRAL}/astral*.jar ]
     do
       read -p "ASTRAL executable file is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/Astral-5.7.8):      " DIR_ASTRAL_TEMP
       DIR_ASTRAL=$(realpath $(echo ${DIR_ASTRAL_TEMP} | sed "s/'//g"))
     done
  echo "ASTRAL ...... OK"

  #check IQ-TREE
  if [ $(which iqtree2) ]
    then
      echo "IQ-TREE ...... OK"
      EXE_IQTREE=$(which iqtree2)
      DIR_IQTREE=${EXE_IQTREE%/*}
    else
      until [ -x $EXE_IQTREE ]
        do
          read -p "IQ-TREE is not found. Please input the location of its executive file (e.g. /usr/local/bin/iqtree2):      " EXE_IQTREE_TEMP
          EXE_IQTREE=$(realpath $(echo ${EXE_IQTREE_TEMP} | sed "s/'//g"))
          DIR_IQTREE=${EXE_IQTREE%/*}
        done
      echo "IQ-TREE ...... OK"
  fi

  #input supermatrix-based (T1) phylogeny
  echo
  read -p "Please input the Newick-formated topology (T1) from the supermatrix-based phylogenies:    " T1_TEMP
  T1=$(realpath $(echo ${T1_TEMP} | sed "s/'//g"))
    until [ -s ${T1} ]
     do
       read -p "Please input the Newick-formated topology (T1) from the supermatrix-based phylogenies:    " T1_TEMP
       T1=$(realpath $(echo ${T1_TEMP} | sed "s/'//g"))
     done
  echo "topology T1 ...... OK"

  #input supertree-based (T2) phylogeny
  echo
  read -p "Please input the Newick-formated topology (T2) from the supertree-based phylogenies:    " T2_TEMP
  T2=$(realpath $(echo ${T2_TEMP} | sed "s/'//g"))
      until [ -s ${T2} ]
     do
       read -p "Please input the Newick-formated topology (T2) from the supermatrix-based phylogenies:    " T2_TEMP
       T1=$(realpath $(echo ${T1_TEMP} | sed "s/'//g"))
     done
  echo "topology T2 ...... OK"

  mkdir -p ${DIR_OUTPUT}/temp/GLS ${DIR_OUTPUT}/temp/GQS ${DIR_OUTPUT}/consistent_gene_results

  #Estimate site-wise loglikelihood for alternative hypotheses

  #construct a supermatrix
  echo
  echo "construct a supermatrix ......"
  cd ${DIR_INPUT_ALIGN}
  ${DIR_PHYKIT}/phykit cc -a ${DIR_OUTPUT}/loci.list -p GLS
  sed -i 's/.*,/LG,/g' GLS.partition   #set model in parition file as 'LG' instead of "AUTO"
  mv GLS.* ${DIR_OUTPUT}/temp/GLS && cd ${DIR_OUTPUT}/temp/GLS
  cat ${T1} ${T2} > ML2ASTRAL.tre

  #generate log-likelihoods of a set of trees and site log-likelihoods
  ${DIR_IQTREE}/iqtree2 -s GLS.fa -Q GLS.partition -m LG+F+R4 -z ML2ASTRAL.tre -wsl --threads-max ${THREADS} --prefix ML2ASTRAL
  awk '{print $2}' GLS.partition | sed 's/=/\t/g;s/-/\t/g' > loci.range
  sed -n 2p ML2ASTRAL.sitelh | awk '{for (i=2; i<=NF; i++) print $i}' > T1.sitelh
  sed -n 3p ML2ASTRAL.sitelh | awk '{for (i=2; i<=NF; i++) print $i}' > T2.sitelh

  #site log-likelihoods of each gene upon T1/T2
  echo -e "loci""\t""GLS_T1""\t""GLS_T2""\t""diff_GLS""\t""support" > GLS_table.txt
  for loci in $(cat ${DIR_OUTPUT}/loci.list)
    do
      start=$(awk '/'"${loci}"'/ {print $2}' loci.range)
      end=$(awk '/'"${loci}"'/ {print $3}' loci.range)
      GLS_T1=$(sed -n "${start}, ${end}p" T1.sitelh | awk '{sum+=$1}END{print sum}')
      GLS_T2=$(sed -n "${start}, ${end}p" T2.sitelh | awk '{sum+=$1}END{print sum}')
      diff_GLS=$(echo | awk '{print '"${GLS_T1}"' - '"${GLS_T2}"'}')

      #determine the topology hypothesis supported for each locus
      if [ $(echo "${diff_GLS} > 0" | bc) -eq 1 ]; then
        support="T1"
      elif [ $(echo "${diff_GLS} < 0" | bc) -eq 1 ]; then
        support="T2"
      else
        support="none"
      fi

      #construct GLS table
      echo -e "${loci}""\t""${GLS_T1}""\t""${GLS_T2}""\t""${diff_GLS}""\t""${support}" >> GLS_table.txt
    done

  mv GLS_table.txt ML2ASTRAL.sitelh loci.range ${DIR_OUTPUT}/consistent_gene_results/


  #estimate the quartet scores for alternative hypotheses
  cd ${DIR_OUTPUT}/temp/GQS
  ${DIR_PHYKIT}/phykit tl ${T1} > tips #tip labels of the species tree

  export DIR_PHYKIT DIR_ASTRAL DIR_INPUT_TREE T1 T2
  temp_fun() {
    #generate the species tree of the identical tip lables with the gene tree
    ${DIR_PHYKIT}/phykit tl ${DIR_INPUT_TREE}/"$1"* > tips."$1"
    grep -v -f tips."$1" tips > tips.missing."$1"
    ${DIR_PHYKIT}/phykit prune ${T1} tips.missing."$1" -o T1.prune.tre."$1"
    ${DIR_PHYKIT}/phykit prune ${T2} tips.missing."$1" -o T2.prune.tre."$1"

    #calculate quartet score for alternative hypotheses
    java -jar ${DIR_ASTRAL}/astral*.jar -i ${DIR_INPUT_TREE}/"$1"* -q T1.prune.tre."$1" 2> log.txt.T1."$1"
    GQS_T1=$(awk '/^Final quartet score/ {print $5}' log.txt.T1."$1")
    java -jar ${DIR_ASTRAL}/astral*.jar -i ${DIR_INPUT_TREE}/"$1"* -q T2.prune.tre."$1" 2> log.txt.T2."$1"
    GQS_T2=$(awk '/^Final quartet score/ {print $5}' log.txt.T2."$1")
    diff_GQS=$(echo | awk '{print '"${GQS_T1}"' - '"${GQS_T2}"'}')

    #determine the topology hypothesis supported for each locus
    if [ $(echo "${diff_GQS} > 0" | bc) -eq 1 ]; then
      support="T1"
    elif [ $(echo "${diff_GQS} < 0" | bc) -eq 1 ]; then
      support="T2"
    else
      support="none"
    fi

    echo -e "$1""\t""${GQS_T1}""\t""${GQS_T2}""\t""${diff_GQS}""\t""${support}" >> GQS."$1"
  }
  export -f temp_fun
  cat ${DIR_OUTPUT}/loci.list | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 temp_fun %

  #construct GQS table
  cat GQS.* > GQS_table.txt
  sed -i "1i "loci"\t"GQS_T1"\t"GQS_T2"\t"diff_GQS"\t"support"" GQS_table.txt
  mv GQS_table.txt ${DIR_OUTPUT}/consistent_gene_results/

  
  #get the consistent loci
  cd ${DIR_OUTPUT}/consistent_gene_results/
  for loci in $(cat ${DIR_OUTPUT}/loci.list)
    do
      GLS=$(awk '/^'"${loci}"'/ {print $5}' GLS_table.txt)
      GQS=$(awk '/^'"${loci}"'/ {print $5}' GQS_table.txt)
      [[ ${GLS} == ${GQS} && ${GLS} != "none" ]] && (echo ${loci} >> ${DIR_OUTPUT}/loci.remaining)
    done
  rm -rf ${DIR_OUTPUT}/temp

fi


#cpoy corresponding alignments and gene trees
if [ ${FILTER_METHOD} -ne 7 ]; then
  cd ${DIR_OUTPUT}
  for loci in $(cat loci.remaining); do cp ${DIR_INPUT_ALIGN}/${loci} alignments_remaining/; cp ${DIR_INPUT_TREE}/${loci}* gene_trees_remaining/; done

  #generate basic statistics
  cd ${DIR_OUTPUT}/alignments_remaining
  ${DIR_PHYKIT}/phykit cc -a ../loci.remaining -p matrix #generate matrix  
  LOCI_TOTAL=$(cat ../loci.list | wc -l) #number of loci
  LOCI_REMAINING=$(cat ../loci.remaining | wc -l) #number of remaining loci
  SITES_REMAINING=$(tail -n 1 matrix.partition | cut -d "-" -f2) #total sites of remaining loci
  AVERAGE=$(echo "scale=2;(${SITES_REMAINING}/${LOCI_REMAINING})"|bc) #average length of remaining loci
  rm matrix* && cd ${DIR_OUTPUT}

  echo
  echo "Among ${LOCI_TOTAL} loci, ${LOCI_REMAINING} of them and ${SITES_REMAINING} sites have been preserved. ${LOCI_REMAINING} loci have an average length of $AVERAGE." | tee -a ${DIR_OUTPUT}/summary.statistics
  echo
  echo "The file loci.remaining contains all the loci list of remaining alignments."
  echo
  echo "All the remaining alignments and gene trees are deposited in the folder 'alignments_remaining' and 'gene_trees_remaining', respectively."

fi





