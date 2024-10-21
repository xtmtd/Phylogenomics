#!/bin/bash
#2021.08.04 wrote by ZF
#2023.06.29 revised by DSY
#2023.10.15 revised by DSY
#2024.07.09 revised by DSY

#Trim alignments using trimal, BMGE, ClipKIT and PhyKIT
#Type "bash trimming_alignments.sh"
#parallel and at least one of the three trimming tools should be installed
#ClipKIT were installed by contructing CONDA virtual environment 'clipkit' and activated by 'source activate clipkit'


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


#Check the threads can be used
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
  until [ ${THREADS} -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done


#Check the type of input alignments
read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide   " ALIGN_TYPE
  until [ ${ALIGN_TYPE} -gt 0 -a ${ALIGN_TYPE} -lt 3 ]
    do
      read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide      " ALIGN_TYPE
    done


#Check the trimming tool can be used
read -p "Please input the option for trimming tool: 1. trimAl; 2. BMGE; 3. ClipKIT   " TRIM_METHOD
  until [ ${TRIM_METHOD} -gt 0 -a ${TRIM_METHOD} -lt 4 ]
    do
      read -p "Please input the option for trimming tool: 1. trimAl; 2. BMGE; 3. ClipKIT      " TRIM_METHOD
    done


#check the trimming tools
#Check the trimal installation and the method used by trimal
if [ ${TRIM_METHOD} == "1" ]; then
  #check trimAL installation
  if [ $(which trimal) ]
    then
      echo "trimal ...... OK"
      EXE_TRIMAL=$(which trimal)
      DIR_TRIMAL=${EXE_TRIMAL%/*}
    else
      until [ -x ${DIR_TRIMAL}/trimal ]
        do
          read -p "trimal is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/trimal-1.4.1/source):      " DIR_TRIMAL_TEMP
          DIR_TRIMAL=$(realpath $(echo ${DIR_TRIMAL_TEMP} | sed "s/'//g"))
        done
      echo "trimal ...... OK"
  fi

  read -p "Please input the option for trimming strategy: 1. -automated1; 2. -gappyout   " TRIMAL_METHOD
    until [ ${TRIMAL_METHOD} -gt 0 -a ${TRIMAL_METHOD} -lt 3 ]
      do
        read -p "Please input the option for trimming strategy: 1. -automated1; 2. -gappyout      " TRIMAL_METHOD
      done

#Check the BMGE installation and the method used by BMGE
elif [ ${TRIM_METHOD} == "2" ]; then
  #check BMGE installation
  if [ $(which BMGE.jar) ]
    then
      echo "BMGE ...... OK"
      EXE_BMGE=$(which BMGE.jar)
      DIR_BMGE=${EXE_BMGE%/*}
    else
      until [ -s ${DIR_BMGE}/BMGE.jar ]
        do
          read -p "BMGE is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/BMGE-1.12):      " DIR_BMGE_TEMP
          DIR_BMGE=$(realpath $(echo ${DIR_BMGE_TEMP} | sed "s/'//g"))
        done
      echo "BMGE ...... OK"
  fi

  read -p "Please input the option for trimming strategy: 1. default; 2. stringent (-h 0.4)   " BMGE_METHOD
    until [ ${BMGE_METHOD} -gt 0 -a ${BMGE_METHOD} -lt 3 ]
      do
        read -p "Please input the option for trimming strategy: 1. default; 2. stringent (-h 0.4)      " BMGE_METHOD
      done

#Check the ClipKIT installation and the method used by ClipKIT
elif [ ${TRIM_METHOD} == "3" ]; then
  #check ClipKIT installation
  if [ $(which clipkit) ]
    then
      echo "clipkit ...... OK"
      EXE_CLIPKIT=$(which clipkit)
      DIR_CLIPKIT=${EXE_CLIPKIT%/*}
    else
      until [ -x ${DIR_CLIPKIT}/clipkit ]
        do
          read -p "ClipKIT is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/ClipKIT/):      " DIR_CLIPKIT_TEMP
          DIR_CLIPKIT=$(realpath $(echo ${DIR_CLIPKIT_TEMP} | sed "s/'//g"))
        done
      echo "ClipKIT ...... OK"
  fi

  #check the method used by ClipKIT
  read -p "Please input the option for trimming strategy: 1. smart-gap (dynamic determination of gaps threshold); 2. kpi-smart-gap (a combination of keeping parsimony informative sites and smart-gap-based trimming); 3. kpic-smart-gap (a combination of keeping parismony informative and constant sites and smart-gap-based trimming   " CLIPKIT_METHOD
    until [ ${CLIPKIT_METHOD} -gt 0 -a ${CLIPKIT_METHOD} -lt 4 ]
      do
        read -p "Please input the option for trimming strategy: 1. smart-gap (dynamic determination of gaps threshold); 2. kpi-smart-gap (a combination of keeping parsimony informative sites and smart-gap-based trimming); 3. kpic-smart-gap (a combination of keeping parismony informative and constant sites and smart-gap-based trimming      " CLIPKIT_METHOD
      done

fi


#input the name of input directory
read -p "Please input the folder names (with its path) containing amino acid/nucleotide alignments, for example '3-align/faa':      " INPUT1_TEMP
INPUT=$(realpath $(echo ${INPUT1_TEMP} | sed "s/'//g"))


#input the name of output directory
read -p "Please input the name of output directory:      " DIR_OUTPUT_TEMP
DIR_OUTPUT=$(echo ${DIR_OUTPUT_TEMP} | sed "s/'//g")
test -d ${DIR_OUTPUT} || mkdir ${DIR_OUTPUT}
cd ${DIR_OUTPUT}


#trimal
if [ ${TRIM_METHOD} == "1" ]; then
  mkdir -p ${DIR_OUTPUT}/trimal

  #trimming
  ls ${INPUT}/ | cut -d "." -f1 > ${DIR_OUTPUT}/loci.list

  if [ ${TRIMAL_METHOD} == "1" ]; then
    cat ${DIR_OUTPUT}/loci.list | ${DIR_PARALLEL}/parallel -j ${THREADS} ${DIR_TRIMAL}/trimal -in ${INPUT}/{}.* -out ${DIR_OUTPUT}/trimal/{}.fas -automated1
  elif [ ${TRIMAL_METHOD} == "2" ]; then
    cat ${DIR_OUTPUT}/loci.list | ${DIR_PARALLEL}/parallel -j ${THREADS} ${DIR_TRIMAL}/trimal -in ${INPUT}/{}.* -out ${DIR_OUTPUT}/trimal/{}.fas -gappyout
  fi

  #trim correspoding nucleotide sequences
  if [ ${ALIGN_TYPE} == "1" ]; then
    echo
    read -p "Do you want to trim the the corresponding nucleotide sequences? 1. Yes; 2. No:      " NUC_TRIM
    until [ ${NUC_TRIM} -gt 0 -a ${NUC_TRIM} -lt 3 ]
      do
        read -p "Do you want to trim the the corresponding nucleotide sequences? 1. Yes; 2. No:      " NUC_TRIM
      done
  fi

  if [[ ${NUC_TRIM} == "1" ]]; then
    #input the name of input nucleotide directory
    read -p "Please input the folder names (with its path) containing unaligned nucleotide alignments, for example '2-loci_filter/fna':      " INPUT2_TEMP
    INPUT_NUC=$(realpath $(echo ${INPUT2_TEMP} | sed "s/'//g"))
    mv ${DIR_OUTPUT}/trimal ${DIR_OUTPUT}/trimal_faa
    mkdir -p ${DIR_OUTPUT}/trimal_fna

    #trimming nucleotide alignments
    if [ ${TRIMAL_METHOD} == "1" ]; then
      cat ${DIR_OUTPUT}/loci.list | ${DIR_PARALLEL}/parallel -j ${THREADS} ${DIR_TRIMAL}/trimal -in ${INPUT}/{}.* -out ${DIR_OUTPUT}/trimal_fna/{}.fas -automated1 -backtrans ${INPUT_NUC}/{}.*
    elif [ ${TRIMAL_METHOD} == "2" ]; then
      cat ${DIR_OUTPUT}/loci.list | ${DIR_PARALLEL}/parallel -j ${THREADS} ${DIR_TRIMAL}/trimal -in ${INPUT}/{}.* -out ${DIR_OUTPUT}/trimal_fna/{}.fas -gappyout -backtrans ${INPUT_NUC}/{}.*
    fi

  else
    echo
  fi

  #detect null alignments
  cd ${DIR_OUTPUT}
  for loci in $(cat loci.list)
    do 
      a=$(sed -n 2p $(ls trimal*/${loci}.* | head -n 1)) 
      b=$(echo ${a:0:1})
      if [ "${b}" == ">" ]; then
       echo ${loci} >> missing.loci
      fi
    done

  #delete null alignments
  if [ -f missing.loci ]; then
    for missing in $(cat missing.loci); do rm trimal*/${missing}.fas; done
    ls trimal*/ | grep -v "trimal" | awk NF | sort -u > loci.list
  fi

fi



#BMGE
if [ ${TRIM_METHOD} == "2" ]; then
  test ${ALIGN_TYPE} == "1" && mkdir -p ${DIR_OUTPUT}/bmge_faa || mkdir -p ${DIR_OUTPUT}/bmge_fna

  #trimming
  ls ${INPUT}/ > ${DIR_OUTPUT}/loci.list

  if [[ ${BMGE_METHOD} == "1" && ${ALIGN_TYPE} == "1" ]]; then
    cat ${DIR_OUTPUT}/loci.list | ${DIR_PARALLEL}/parallel -j ${THREADS} java -jar ${DIR_BMGE}/BMGE.jar -i ${INPUT}/{} -t AA -of ${DIR_OUTPUT}/bmge_faa/{}   
  elif [[ ${BMGE_METHOD} == "2" && ${ALIGN_TYPE} == "1" ]]; then
    cat ${DIR_OUTPUT}/loci.list | ${DIR_PARALLEL}/parallel -j ${THREADS} java -jar ${DIR_BMGE}/BMGE.jar -i ${INPUT}/{} -t AA -m BLOSUM90 -h 0.4 -of ${DIR_OUTPUT}/bmge_faa/{}
  elif [[ ${BMGE_METHOD} == "1" && ${ALIGN_TYPE} == "2" ]]; then
    cat ${DIR_OUTPUT}/loci.list | ${DIR_PARALLEL}/parallel -j ${THREADS} java -jar ${DIR_BMGE}/BMGE.jar -i ${INPUT}/{} -t DNA -of ${DIR_OUTPUT}/bmge_fna/{}
  elif [[ ${BMGE_METHOD} == "2" && ${ALIGN_TYPE} == "2" ]]; then
    cat ${DIR_OUTPUT}/loci.list | ${DIR_PARALLEL}/parallel -j ${THREADS} java -jar ${DIR_BMGE}/BMGE.jar -i ${INPUT}/{} -t DNA -m DNAPAM50 -h 0.4 -of ${DIR_OUTPUT}/bmge_fna/{}
  fi

  #detect null alignments
  cd ${DIR_OUTPUT}
  for loci in $(cat loci.list)
    do 
      a=$(sed -n 2p $(ls bmge*/${loci} | head -n 1)) 
      b=$(echo ${a:0:1})
      if [ "${b}" == ">" ]; then
       echo ${loci} >> missing.loci
      fi
    done

  #delete null alignments
  if [ -f missing.loci ]; then
    for missing in $(cat missing.loci); do rm bmge*/${missing}; done
    ls bmge*/ | grep -v "bmge" | awk NF | sort -u > loci.list
  fi

fi


#clipkit
if [ ${TRIM_METHOD} == "3" ]; then
  mkdir -p ${DIR_OUTPUT}/clipkit

  #trimming
  ls ${INPUT}/ | cut -d "." -f1 > loci.list

  if [ ${CLIPKIT_METHOD} == "1" ]; then
    cat loci.list | ${DIR_PARALLEL}/parallel -j ${THREADS} ${DIR_CLIPKIT}/clipkit ${INPUT}/{}.* -o ${DIR_OUTPUT}/clipkit/{}.fas -l -m smart-gap
  elif [ ${CLIPKIT_METHOD} == "2" ]; then
    cat loci.list | ${DIR_PARALLEL}/parallel -j ${THREADS} ${DIR_CLIPKIT}/clipkit ${INPUT}/{}.* -o ${DIR_OUTPUT}/clipkit/{}.fas -l -m kpi-smart-gap
  elif [ ${CLIPKIT_METHOD} == "3" ]; then
    cat loci.list | ${DIR_PARALLEL}/parallel -j ${THREADS} ${DIR_CLIPKIT}/clipkit ${INPUT}/{}.* -o ${DIR_OUTPUT}/clipkit/{}.fas -l -m kpic-smart-gap
  fi

  #trim correspoding nucleotide sequences
  if [ ${ALIGN_TYPE} == "1" ]; then
    echo
    read -p "Do you want to trim the the corresponding nucleotide sequences? 1. Yes; 2. No:      " NUC_TRIM
    until [ ${NUC_TRIM} -gt 0 -a ${NUC_TRIM} -lt 3 ]
      do
        read -p "Do you want to trim the the corresponding nucleotide sequences? 1. Yes; 2. No:      " NUC_TRIM
      done
  fi

  if [[ ${NUC_TRIM} == "1" ]]; then
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
	
    #input the name of input nucleotide directory
    read -p "Please input the folder names (with its path) containing nucleotide alignments, for example '3-faa_align/fna':      " INPUT2_TEMP
    INPUT_NUC=$(realpath $(echo ${INPUT2_TEMP} | sed "s/'//g"))
    mv ${DIR_OUTPUT}/clipkit ${DIR_OUTPUT}/clipkit_faa
    mkdir -p ${DIR_OUTPUT}/clipkit_fna

    #trimming nucleotide alignments
    cd ${DIR_OUTPUT}
    export DIR_PHYKIT INPUT_NUC
    TEMP_fun() {
      ${DIR_PHYKIT}/phykit p2n -p clipkit_faa/"$1".fas -n ${INPUT_NUC}/"$1".* -c clipkit_faa/"$1".*.log > clipkit_fna/"$1".fas
    }
    export -f TEMP_fun
    cat ../loci.list | ${DIR_PARALLEL}/parallel -j ${THREADS} TEMP_fun
    rm clipkit_faa/*.log

  #detect null alignments
  for loci in $(cat ../loci.list)
    do 
      a=$(sed -n 2p $(ls clipkit*/${loci}.* | head -n 1)) 
      b=$(echo ${a:0:1})
      if [ "${b}" == ">" ]; then
       echo ${loci} >> missing.loci
      fi
    done

  #delete null alignments
  if [ -f missing.loci ]; then
    for missing in $(cat missing.loci); do rm clipkit*/${missing}.fas; done
    ls clipkit*/ | grep -v "clipkit" | awk NF | sort -u > loci.list
  fi

  else
    rm ${DIR_OUTPUT}/clipkit/*.log

  #detect null alignments
  for loci in $(cat loci.list)
    do 
      a=$(sed -n 2p $(ls ${DIR_OUTPUT}/clipkit*/${loci}.* | head -n 1)) 
      b=$(echo ${a:0:1})
      if [ "${b}" == ">" ]; then
       echo ${loci} >> missing.loci
      fi
    done

  #delete null alignments
  if [ -f missing.loci ]; then
    for missing in $(cat missing.loci); do rm clipkit*/${missing}.fas; done
    ls ${DIR_OUTPUT}/clipkit*/ | grep -v "clipkit" | awk NF | sort -u > loci.list
  fi

  fi

fi





