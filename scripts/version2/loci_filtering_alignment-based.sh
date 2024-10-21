#!/bin/bash
#2021.08.04 by ZF v1
#2023.01.13 by DSY v2
#2023.07.17 by DSY v3
#2023.08.28 by DSY v4
#2023.10.22 by ZF

#filter loci by detecting alignment length, number of parsimony-informative sites, percentage of parsimony-informative sites in the alignment, GC content, compositional heterogeneity (RCV, Relative composition variability; nRCFV, normalised Relative Compositional Frequency Variation), average pairwise identity, likelihood mapping criteria and TAPER.
#Type 'bash loci_filtering_alignment-based.sh'
#parallel, PhyKIT, seqkit, csvtk, TAPER, Julia, RCFV_Reader, FAST and IQ-TREE may be required


#check Phykit
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


#check Seqkit
if [ $(which seqkit) ]
    then
      echo "SEQKIT ...... OK"
      EXE_SEQKIT=$(which seqkit)
      DIR_SEQKIT=${EXE_SEQKIT%/*}
    else
      until [ -x ${DIR_SEQKIT}/seqkit ]
  do
    read -p "SEQKIT is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_SEQKIT_TEMP
    DIR_SEQKIT=$(realpath $(echo ${DIR_SEQKIT_TEMP} | sed "s/'//g"))
  done
      echo "SEQKIT ...... OK"
fi


#check csvtk
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
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
  until [ ${THREADS} -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done


DIR_CURR=$(echo ${PWD})

#input the name of input directory
read -p "Please input the name of input directory containing all alignments, e.g. 4-trim/clipkit-kpic:      " DIR_INPUT_TEMP
DIR_INPUT=$(realpath $(echo ${DIR_INPUT_TEMP} | sed "s/'//g"))


#input the name of output directory
read -p "Please input the name of output directory, or an existing directory:      " DIR_OUTPUT_TEMP
DIR_OUTPUT_TEMP1=$(echo ${DIR_OUTPUT_TEMP} | sed "s/'//g")
test -d ${DIR_OUTPUT_TEMP1} || mkdir -p ${DIR_OUTPUT_TEMP1}
cd ${DIR_OUTPUT_TEMP1} && DIR_OUTPUT=$(echo ${PWD}) && cd ${DIR_CURR}


#Check the loci filtering method can be used
read -p "Please input the option for alignment-based strategy for loci filtering: 1. alignment length; 2. number of parsimony-informative sites; 3. percentage of parsimony-informative sites in the alignment; 4. GC content; 5. RCFV (Relative Compositional Frequency Variation); 6. nRCFV (normalised Relative Compositional Frequency Variation); 7. average pairwise identity; 8. likelihood mapping; 9. symmetry tests against SRH hypotheses; 10. TAPER     " FILTER_METHOD
  until [[ ${FILTER_METHOD} -gt 0 && ${FILTER_METHOD} -lt 11 ]]
    do
      read -p "Please input the option for alignment-based strategy for loci filtering: 1. alignment length; 2. number of parsimony-informative sites; 3. percentage of parsimony-informative sites in the alignment; 4. GC content; 5. RCFV (Relative Compositional Frequency Variation); 6. nRCFV (normalised Relative Compositional Frequency Variation); 7. average pairwise identity; 8. likelihood mapping; 9. symmetry tests against SRH hypotheses; 10. TAPER     " FILTER_METHOD
    done


#generate output folders and loci list
mkdir -p ${DIR_OUTPUT}/alignments_remaining ${DIR_OUTPUT}/temp
cd ${DIR_INPUT} && ls * > ${DIR_OUTPUT}/loci.list #generate loci list
cd ${DIR_OUTPUT}


#site length-based filtering
if [ ${FILTER_METHOD} -lt 4 ]; then
  echo
  echo "Calculating alignment length statistics ......"

  #calculate basic statistics for alignments
  export DIR_INPUT DIR_PHYKIT DIR_OUTPUT
  temp_fun() {
    echo -e "$1""\t""$(${DIR_PHYKIT}/phykit pis ${DIR_INPUT}/$1)" > ./temp/$1
  }
  export -f temp_fun
  cat loci.list | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 temp_fun %
  cat temp/* | sort -k3n -k2n > length_site.txt
  rm -rf temp

  #generate the frequency distribution plot and help to determine the threshold
  if [ ${FILTER_METHOD} -eq 1 ]; then
    awk '{print $3}' length_site.txt | sort -k1n | awk -v OFS="\t" '{print NR, $0}' | sed '1i\Loci_number\tNumber of sites' | csvtk -t plot line -x Loci_number -y "Number of sites" --format pdf > Alignment_length.pdf
  elif [ ${FILTER_METHOD} -eq 2 ]; then
    awk '{print $2}' length_site.txt | sort -k1n | awk -v OFS="\t" '{print NR, $0}' | sed '1i\Loci_number\tNumber of parsimony-informative sites' | csvtk -t plot line -x Loci_number -y "Number of parsimony-informative sites" --format pdf > Parsimony_length.pdf
  elif [ ${FILTER_METHOD} -eq 3 ]; then
    awk '{print $4}' length_site.txt | sort -k1n | awk -v OFS="\t" '{print NR, $0}' | sed '1i\Loci_number\tPercentage of parsimony-informative sites' | csvtk -t plot line -x Loci_number -y "Percentage of parsimony-informative sites" --format pdf > Percentage_of_parsimony_sites.pdf
  else
    echo
  fi

  #input the threshold 
  echo
  echo "Read the file length_site.txt in the output folder to determine the threshold. The 2nd~4th columns represent the number of parsimony-informative sites, the alignment length, the percentage of parsimony-informative sites in the alignment. The distribution plot (XXX.pdf) may also help to determine the cutoff values."
  echo
  read -p "Please input the threshold (the minimum value):      " THRESHOLD

  #generate the list of remaining loci
  if [ ${FILTER_METHOD} -eq 1 ]; then
    awk '$3>='"${THRESHOLD}"' {print $1}' length_site.txt > loci.remaining
  elif [ ${FILTER_METHOD} -eq 2 ]; then
    awk '$2>='"${THRESHOLD}"' {print $1}' length_site.txt > loci.remaining
  elif [ ${FILTER_METHOD} -eq 3 ]; then
    awk '$4>='"${THRESHOLD}"' {print $1}' length_site.txt > loci.remaining
  else
    echo
  fi


#GC content
elif [ ${FILTER_METHOD} -eq 4 ]; then  
  echo

  #calculate GC content
  echo "Calculating GC content ......"

  export DIR_INPUT DIR_PHYKIT DIR_OUTPUT
  temp_fun() {
    echo -e "$1""\t""$(${DIR_PHYKIT}/phykit gc ${DIR_INPUT}/$1)" > ./temp/$1
  }
  export -f temp_fun
  cat loci.list | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 temp_fun %
  cat temp*/* | sort -k2n > GC.txt
  rm -rf temp

  #generate the frequency distribution plot and help to determine the threshold
  awk '{print $2}' GC.txt | sort -k1n | awk -v OFS="\t" '{print NR, $0}' | sed '1i\Loci_number\tGC content' | csvtk -t plot line -x Loci_number -y "GC content" --format pdf > GC.pdf

  #input the threshold 
  echo
  echo "Read the file GC.txt in the output folder to determine the GC percentage threshold. Very high (e.g. > 0.8) and/or very low (often < 0.2) GC values are thought to be unreliable. The distribution plot GC.pdf may also help to determine the cutoff values."
  echo
  read -p "Please input the minimum GC percentage threshold:      " MIN_GC
  echo
  read -p "Please input the maximum GC percentage threshold:      " MAX_GC
  
  #generate the list of remaining loci
  awk '$2>='"${MIN_GC}"'&&$2<='"${MAX_GC}"' {print $1}' GC.txt > loci.remaining


#RCFV/nRCFV
elif [ ${FILTER_METHOD} -eq 5 -o ${FILTER_METHOD} -eq 6 ]; then  
#check RCFVReader_v1.pl
if [ $(which RCFVReader_v1.pl) ]
    then
      echo "RCFVReader_v1.pl ...... OK"
      EXE_RCFVREADER=$(which RCFVReader_v1.pl)
      DIR_RCFVREADER=${EXE_RCFVREADER%/*}
    else
      until [ -x ${DIR_RCFVREADER}/RCFVReader_v1.pl ]
        do
          read -p "RCFVReader_v1.pl is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/RCFV_Reader):      " DIR_RCFVREADER_TEMP
          DIR_RCFVREADER=$(realpath $(echo ${DIR_RCFVREADER_TEMP} | sed "s/'//g"))
        done
      echo "RCFVReader_v1.pl ...... OK"
fi

  echo
  #Check the type of input alignments
  read -p "Please input the sequence type of input alignments: 1. amino acid; 2. nucleotide   " ALIGN_TYPE
    until [[ ${ALIGN_TYPE} -gt 0 && ${ALIGN_TYPE} -lt 3 ]]
      do
        read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide      " ALIGN_TYPE
      done

  #calculate RCFV
  echo "Calculating RCFV/nRCFV ......"

  cd ${DIR_OUTPUT}/temp/
  export DIR_INPUT DIR_OUTPUT ALIGN_TYPE DIR_RCFVREADER
  temp_fun() {
    test ${ALIGN_TYPE} -eq 1 && perl ${DIR_RCFVREADER}/RCFVReader*.pl protein ${DIR_INPUT}/$1 $1 || perl ${DIR_RCFVREADER}/RCFVReader*.pl dna ${DIR_INPUT}/$1 $1
    RCFV=$(awk '/^RCFV/ {print $2}' "$1".RCFV.txt)
    nRCFV=$(awk '/^nRCFV/ {print $2}' "$1".RCFV.txt)
    echo -e "$1""\t""${RCFV}""\t""${nRCFV}" > $1
    rm "$1"*txt
  }
  export -f temp_fun
  cat ../loci.list | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 temp_fun %
  cd ..
  cat temp/* | sort -k2n -k3n > RCFV.txt
  rm -rf temp

  #generate the frequency distribution plot and help to determine the threshold
  if [ ${FILTER_METHOD} -eq 5 ]; then
    awk '{print $2}' RCFV.txt | sort -k1n | awk -v OFS="\t" '{print NR, $0}' | sed '1i\Loci_number\tRCFV' | csvtk -t plot line -x Loci_number -y "RCFV" --format pdf > RCFV.pdf
  elif [ ${FILTER_METHOD} -eq 6 ]; then
    awk '{print $3}' RCFV.txt | sort -k1n | awk -v OFS="\t" '{print NR, $0}' | sed '1i\Loci_number\tnRCFV' | csvtk -t plot line -x Loci_number -y "nRCFV" --format pdf > nRCFV.pdf
  else
    echo
  fi

  #input the threshold 
  echo
  echo "Read the file RCFV.txt in the output folder to determine the RCFV/nRCFV threshold. Lower values are thought to be desirable because they represent a lower composition bias in an alignment. The distribution plot (XXX.pdf) may also help to determine the cutoff values."
  echo
  read -p "Please input the threshold (the maximum value):      " THRESHOLD

  #generate the list of remaining loci
  if [ ${FILTER_METHOD} -eq 5 ]; then
    awk '$2<='"${THRESHOLD}"' {print $1}' RCFV.txt > loci.remaining
  elif [ ${FILTER_METHOD} -eq 6 ]; then
    awk '$3<='"${THRESHOLD}"' {print $1}' RCFV.txt > loci.remaining
  else
    echo
  fi


#average pairwise identity
elif [ ${FILTER_METHOD} -eq 7 ]; then  
  echo

  #calculate average pairwise identity
  echo "Calculating average pairwise identity ......"

  export DIR_INPUT DIR_PHYKIT DIR_OUTPUT
  temp_fun() {
    echo -e "$1""\t""$(${DIR_PHYKIT}/phykit pi ${DIR_INPUT}/$1 | awk 'NR == 1 {print $2}')" > ./temp/$1
  }
  export -f temp_fun
  cat loci.list | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 temp_fun %
  cat temp/* | sort -k2n > API.txt
  rm -rf temp

  #generate the frequency distribution plot and help to determine the threshold
  awk '{print $2}' API.txt | sort -k1n | awk -v OFS="\t" '{print NR, $0}' | sed '1i\Loci_number\taverage pairwise identity' | csvtk -t plot line -x Loci_number -y "average pairwise identity" --format pdf > API.pdf

  #input the threshold 
  echo
  echo "Read the file API.txt in the output folder to determine the rate threshold. Very low values indicate that those 'fast' genes may be more suitable for shallow phylogenies. The distribution plot API.pdf may also help to determine the cutoff values."
  echo
  read -p "Please input the minimum API threshold:      " MIN_API
  echo
  read -p "Please input the maximum API threshold:      " MAX_API
  
  #generate the list of remaining loci
  awk '$2>='"${MIN_API}"'&&$2<='"${MAX_API}"' {print $1}' API.txt > loci.remaining


#likelihood mapping
elif [ ${FILTER_METHOD} -eq 8 ]; then  
  #check IQ-TREE
  if [ $(which iqtree2) ]
    then
      echo "IQ-TREE ...... OK"
      EXE_IQTREE=$(which iqtree2)
      DIR_IQTREE=${EXE_IQTREE%/*}
    else
      until [ -x ${EXE_IQTREE} ]
        do
          read -p "IQ-TREE is not found. Please input the location of its executive file (e.g. /usr/local/bin/iqtree2):      " EXE_IQTREE_TEMP
          EXE_IQTREE=$(realpath $(echo ${EXE_IQTREE_TEMP} | sed "s/'//g"))
          DIR_IQTREE=${EXE_IQTREE%/*}
        done
      echo "IQ-TREE ...... OK"
  fi

  echo
  #Check the type of input alignments
  read -p "Please input the sequence type of input alignments: 1. amino acid; 2. nucleotide   " ALIGN_TYPE
    until [[ ${ALIGN_TYPE} -gt 0 && ${ALIGN_TYPE} -lt 3 ]]
      do
        read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide      " ALIGN_TYPE
      done


  echo
  #calculate likelihood mapping
  echo "Calculating likelihood mapping ......"
  cp ${DIR_INPUT}/* ${DIR_OUTPUT}/temp/

  if [ ${ALIGN_TYPE} -eq 1 ]; then
  cd ${DIR_OUTPUT}/temp/
  export DIR_INPUT DIR_OUTPUT DIR_IQTREE
  temp_fun() {
    ${DIR_IQTREE}/iqtree2 -s $1 -m LG --lmap 10000 -n 0 -T 1
    echo -e "$1""\t""$(awk '/^Number of fully resolved  quartets/ {print $8/10000}' "$1".iqtree)" > "$1".LMQ
  }
  export -f temp_fun
  cat ../loci.list | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 temp_fun %
  cd ..

  else
  cd ${DIR_OUTPUT}/temp/
  export DIR_INPUT DIR_OUTPUT DIR_IQTREE
  temp_fun() {
    ${DIR_IQTREE}/iqtree2 -s $1 -m GTR --lmap 10000 -n 0 -T 1
    echo -e "$1""\t""$(awk '/^Number of fully resolved  quartets/ {print $8/10000}' "$1".iqtree)" > "$1".LMQ
  }
  export -f temp_fun
  cat ../loci.list | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 temp_fun %
  cd ..

  fi
  cat temp/*.LMQ | sort -k2n > Likelihood_mapping_quartet.txt
  rm -rf temp

  #generate the frequency distribution plot and help to determine the threshold
  awk '{print $2}' Likelihood_mapping_quartet.txt | sort -k1n | awk -v OFS="\t" '{print NR, $0}' | sed '1i\Loci_number\tNumber of fully resolved  quartets' | csvtk -t plot line -x Loci_number -y "Number of fully resolved  quartets" --format pdf > Likelihood_mapping_quartet.pdf

  #input the threshold 
  echo
  echo "Read the file likelihood_mapping_quartet.txt in the output folder to determine the quartet threshold. High values indicate that strong phylogenetic informativeness/signal. The distribution plot 'Likelihood_mapping_quartet.pdf' may also help to determine the cutoff values."
  echo
  read -p "Please input the threshold (the minimum value):      " THRESHOLD

  #generate the list of remaining loci
  awk '$2>='"${THRESHOLD}"' {print $1}' Likelihood_mapping_quartet.txt > loci.remaining


#symmetry tests against SRH hypotheses
#likelihood mapping
elif [ ${FILTER_METHOD} -eq 9 ]; then  
  #check IQ-TREE
  if [ $(which iqtree2) ]
    then
      echo "IQ-TREE ...... OK"
      EXE_IQTREE=$(which iqtree2)
      DIR_IQTREE=${EXE_IQTREE%/*}
    else
      until [ -x ${EXE_IQTREE} ]
        do
          read -p "IQ-TREE is not found. Please input the location of its executive file (e.g. /usr/local/bin/iqtree2):      " EXE_IQTREE_TEMP
          EXE_IQTREE=$(realpath $(echo ${EXE_IQTREE_TEMP} | sed "s/'//g"))
          DIR_IQTREE=${EXE_IQTREE%/*}
        done
      echo "IQ-TREE ...... OK"
  fi

  #generate concatenated matrix
  cd ${DIR_OUTPUT}/temp
  cp -r ${DIR_INPUT}/* .
  ${DIR_PHYKIT}/phykit cc -a ../loci.list -p matrix
  
  #calculate symtest
  echo
  echo "Calculating symmetry tests ......"

  ${DIR_IQTREE}/iqtree2 -s matrix.fa -p matrix.partition --symtest-only -T ${THREADS}
  grep -v "^#" matrix.partition.symtest.csv | cut -d "," -f1,4 | sed "1d;s/,/\t/g" | sort -k2nr > ${DIR_OUTPUT}/symtest.txt
  cd ${DIR_OUTPUT} && rm -rf temp

  #input the threshold 
  echo
  echo "Read the file symtest.txt in the output folder to determine the p-value threshold of symmetry tests, usually 0.01~0.1. Higher values indicate that more loci rejecting SRH hypotheses would be removed."
  echo
  read -p "Please input the threshold:      " THRESHOLD

  #generate the list of remaining loci
  awk '$2>='"${THRESHOLD}"' {print $1}' symtest.txt > loci.remaining


#TAPER check errorous alignment regions
elif [ ${FILTER_METHOD} -eq 10 ]; then  
  #check julia
  if [ $(which julia) ]
    then
      echo "julia ...... OK"
      EXE_JULIA=$(which julia)
      DIR_JULIA=${EXE_JULIA%/*}
    else
      until [ -x ${EXE_JULIA} ]
        do
          read -p "julia is not found. Please input the location of its executive file (e.g. /usr/local/bin/julia):      " EXE_JULIA_TEMP
          EXE_JULIA=$(realpath $(echo ${EXE_JULIA_TEMP} | sed "s/'//g"))
          DIR_JULIA=${EXE_JULIA%/*}
        done
      echo "julia ...... OK"
  fi

  #check TAPER
  until [ -s ${DIR_TAPER}/correction_multi.jl ]
     do
       read -p "correction_multi.jl is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/TAPER-1.0.0):      " DIR_TAPER_TEMP
       DIR_TAPER=$(realpath $(echo ${DIR_TAPER_TEMP} | sed "s/'//g"))
     done
  echo "correction_multi.jl ...... OK"

  #Check the type of input alignments
  read -p "Please input the sequence type of input alignments: 1. amino acid; 2. nucleotide   " ALIGN_TYPE
    until [[ ${ALIGN_TYPE} -gt 0 && ${ALIGN_TYPE} -lt 3 ]]
      do
        read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide      " ALIGN_TYPE
      done

  echo
  #calculate TAPER
  echo "Calculating TAPER ......"

  cd ${DIR_OUTPUT}/alignments_remaining
  export DIR_INPUT DIR_OUTPUT DIR_JULIA DIR_TAPER ALIGN_TYPE
  temp_fun() {
    test ${ALIGN_TYPE} -eq 1 && (${DIR_JULIA}/julia ${DIR_TAPER}/correction_multi.jl ${DIR_INPUT}/$1 > $1) || (${DIR_JULIA}/julia ${DIR_TAPER}/correction_multi.jl -m N -a N ${DIR_INPUT}/$1 > $1)
  }
  export -f temp_fun
  cat ../loci.list | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 temp_fun %
  ls > ../loci.remaining
  cd ${DIR_OUTPUT} && rm -rf temp

fi


#cpoy corresponding alignments
for loci in $(cat loci.remaining); do cp ${DIR_INPUT}/${loci} alignments_remaining/; done

#generate basic statistics
cd ${DIR_OUTPUT}/alignments_remaining
${DIR_PHYKIT}/phykit cc -a ../loci.remaining -p matrix #generate matrix  
LOCI_TOTAL=$(cat ../loci.list | wc -l) #number of loci
LOCI_REMAINING=$(cat ../loci.remaining | wc -l) #number of remaining loci
SITES_REMAINING=$(tail -n 1 matrix.partition | cut -d "-" -f2) #total sites of remaining loci
AVERAGE=$(echo "scale=2;(${SITES_REMAINING}/${LOCI_REMAINING})"|bc) #average length of remaining loci
rm matrix* && cd ${DIR_OUTPUT}

echo
echo "Among ${LOCI_TOTAL} loci, ${LOCI_REMAINING} of them and ${SITES_REMAINING} sites have been preserved. ${LOCI_REMAINING} loci have an average length of ${AVERAGE}." | tee -a ${DIR_OUTPUT}/summary.statistics
echo
echo "The file loci.remaining contains all the file list of remaining alignments."
echo
echo "All the remaining alignments are deposited in the folder 'alignments_remaining'."


