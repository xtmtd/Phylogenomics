#!/bin/bash
#2021.08.04 wrote by ZF, 2023.1.11 modified by DSY
#2023.10.13 by ZF
#2024.05.06 by DSY

#generate the supermaxtrix, partition and occupancy for alignments
#Type 'bash matrix_generation.sh'
#FASconCAT, PhyKIT and Seqkit may be required


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


DIR_CURR=$(echo ${PWD})

#input the name of input directory
read -p "Please input the name of input directory containing all alignments, e.g. 4-trim/clipkit-kpi:      " DIR_INPUT_TEMP
DIR_INPUT=$(realpath $(echo ${DIR_INPUT_TEMP} | sed "s/'//g"))


#input the name of output directory
read -p "Please input the name of output directory, or an existing directory:      " DIR_OUTPUT_TEMP
DIR_OUTPUT_TEMP1=$(echo ${DIR_OUTPUT_TEMP} | sed "s/'//g")
test -d ${DIR_OUTPUT_TEMP1} || mkdir -p ${DIR_OUTPUT_TEMP1}
cd ${DIR_OUTPUT_TEMP1} && DIR_OUTPUT=$(echo ${PWD}) && cd ${DIR_CURR}


#input the prefix for the generated matrix-related files
read -p "Please input the name of a prefix for the generated matrix-related files, e.g. DATASET1:      " PREFIX


#input the taxon occupancy values for the matrices
read -p "Please input the minimum percentage value for taxa occupancy, usually ranging from 50% to 100%, e.g. 50, 75, 90:      " OCCUPANCY
OCCU=$(echo -e ${OCCUPANCY} | awk '{print $1/100}')


#Check the type of input alignments
read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide   " ALIGN_TYPE
  until [ ${ALIGN_TYPE} -gt 0 -a ${ALIGN_TYPE} -lt 3 ]
    do
      read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide      " ALIGN_TYPE
    done


#generate matrix
cd ${DIR_INPUT}
ls > ../${PREFIX}.alignments
${DIR_PHYKIT}/phykit cc -a ../${PREFIX}.alignments -p ${PREFIX}
rm ../${PREFIX}.alignments

mkdir -p ${DIR_OUTPUT}/all ${DIR_OUTPUT}/matrix${OCCUPANCY}/alignments
mv ${PREFIX}.* ${DIR_OUTPUT}/all

#generate the list of loci with taxa occupancy above the threshold
awk '$4>='"${OCCU}"' {print $1}' ${DIR_OUTPUT}/all/*.occupancy > ${DIR_OUTPUT}/matrix${OCCUPANCY}/${PREFIX}.alignments

cd ${DIR_OUTPUT}/matrix${OCCUPANCY}/alignments
for file in $(cat ${DIR_OUTPUT}/matrix${OCCUPANCY}/${PREFIX}.alignments); do cp ${DIR_INPUT}/${file} .; done
$DIR_PHYKIT/phykit cc -a ../${PREFIX}.alignments -p ${PREFIX}
mv  ${PREFIX}.* ${DIR_OUTPUT}/matrix${OCCUPANCY}/
rm -rf ${DIR_OUTPUT}/all


#Whether to recode the alignments
read -p "Whether to recode the alignments: 1. yes; 2. no   " ALIGN_RECOD
until [ ${ALIGN_RECOD} -gt 0 -a ${ALIGN_RECOD} -lt 3 ]
  do
    read -p "Whether to recode the alignments: 1. yes; 2. no   " ALIGN_RECOD
  done

if [ ${ALIGN_RECOD} == "1" ]; then
  #Codes for which recoding scheme to use:
  read -p "Please input a recoding table containing recoded characters, e.g. /home/zf/install/PhyKIT-1.19.4/phykit/recoding_tables/Dayhoff-6.txt:      " RECODING_TEMP
  RECODING=$(realpath $(echo ${RECODING_TEMP} | sed "s/'//g"))
  until [ -s ${RECODING} ]
     do
       read -p "Please input a recoding table containing recoded characters, e.g. /home/zf/install/PhyKIT-1.19.4/phykit/recoding_tables/Dayhoff-6.txt:      " RECODING_TEMP
       RECODING=$(realpath $(echo ${RECODING_TEMP} | sed "s/'//g"))
     done

    phykit recode ${DIR_OUTPUT}/matrix${OCCUPANCY}/${PREFIX}.fa -c ${RECODING} > ${DIR_OUTPUT}/matrix${OCCUPANCY}/${PREFIX}.recoding.fa
  
  else
    echo
  
  fi


#construct matrix with the 3rd sites excluded

#Check the type of input alignments
read -p "For coding sequences of protein-coding genes, one often wants to exlude the 3rd codon positions: 1. exclude 3rd coden sites; 2. keep all   " CODON12
until [ ${CODON12} -gt 0 -a ${CODON12} -lt 3 ]
  do
    read -p "For coding sequences of protein-coding genes, one often wants to exlude the 3rd codon position: 1. exclude 3rd coden sites; 2. keep all   " CODON12
  done

if [ ${CODON12} == "2" ]; then
  echo
else
  #check FASconCAT-G
  if [ $(which FASconCAT-G.pl) ]
    then
      echo "FASconCAT-G.pl ...... OK"
      EXE_FASCONCAT=$(which FASconCAT-G.pl)
      DIR_FASCONCAT=${EXE_FASCONCAT%/*}
    else
      until [ -s $DIR_FASCONCAT/FASconCAT-G*pl ]
        do
          read -p "FASconCAT-G is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/FASconCAT-G-1.05):      " DIR_FASCONCAT_TEMP
          DIR_FASCONCAT=$(realpath $(echo $DIR_FASCONCAT_TEMP | sed "s/'//g"))
        done
      echo "FASconCAT-G.pl ...... OK"
  fi

  mkdir -p ${DIR_OUTPUT}/matrix${OCCUPANCY}/temp
  cd ${DIR_OUTPUT}/matrix${OCCUPANCY}/temp
  cp ${DIR_OUTPUT}/matrix${OCCUPANCY}/alignments/* .

  #generate matrix_cds12
  perl ${DIR_FASCONCAT}/FASconCAT-G*.pl -s -l -d
  mv FcC_supermatrix.fas ../${PREFIX}.cds12.fa
  mv *.txt ../${PREFIX}.cds12.partition
  cd .. && rm -rf temp/
fi

cd ${DIR_OUTPUT}/matrix${OCCUPANCY}

#Move outgroup taxa to the first line
echo
echo "Many phylogenetic tools, such as IQ-TREE, will view the first taxon of the alignment as the 'outgroup' in the final tree file."
echo
read -p "Move the outgroup species to the first one in the alignment matrix? 1. Yes; 2. No:      " MOVE_OUTGROUP
  until [ $MOVE_OUTGROUP -gt 0 -a $MOVE_OUTGROUP -lt 3 ]
    do
      read -p "Move the outgroup species to the first one in the alignment matrix? 1. Yes; 2. No:      " MOVE_OUTGROUP
    done

if [ "$MOVE_OUTGROUP" == "1" ]; then
  #check Seqkit
  if [ $(which seqkit) ]
    then
      echo "SeqKit ...... OK"
      EXE_SEQKIT=$(which seqkit)
      DIR_SEQKIT=${EXE_SEQKIT%/*}
    else
      until [ -x ${DIR_SEQKIT}/seqkit ]
        do
          read -p "SeqKit is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_SEQKIT_TEMP
          DIR_SEQKIT=$(echo ${DIR_SEQKIT_TEMP} | sed "s/'//g")
        done
      echo "SeqKit ...... OK"
  fi

  #input outgrpup name
  read -p "Please input the species name of outgroup taxon, such as 'Zootermopsis_nevadensis':      "  OUTGROUP

  #move outgroup to the first row
  echo ${OUTGROUP} > list.outgroup
  ${DIR_SEQKIT}/seqkit grep -f list.outgroup ${PREFIX}.fa  | ${DIR_SEQKIT}/seqkit seq -u -w 0 > outgroup.fa
  ${DIR_SEQKIT}/seqkit grep -v -f list.outgroup ${PREFIX}.fa | ${DIR_SEQKIT}/seqkit seq -u -w 0 > ingroup.fa
  cat outgroup.fa ingroup.fa > ${PREFIX}.fa

  if [ ${ALIGN_TYPE} == "2" -a ${CODON12} == "1" ]; then
    ${DIR_SEQKIT}/seqkit grep -f list.outgroup ${PREFIX}.cds12.fa  | ${DIR_SEQKIT}/seqkit seq -u -w 0 > outgroup.fa
    ${DIR_SEQKIT}/seqkit grep -v -f list.outgroup ${PREFIX}.cds12.fa | ${DIR_SEQKIT}/seqkit seq -u -w 0 > ingroup.fa
    cat outgroup.fa ingroup.fa > ${PREFIX}.cds12.fa
  else
    echo
  fi
  
  rm outgroup.fa ingroup.fa list.outgroup

else
  echo
fi


#summary thre results
cd ${DIR_OUTPUT}
echo
echo "Individual loci alignments, concatenated matrix and partition file are deposited in the OUTPUT/matrix${OCCUPANCY}."
echo "${OCCUPANCY}% occupancy matrix has $(cat ${DIR_OUTPUT}/matrix${OCCUPANCY}/${PREFIX}.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat ${DIR_OUTPUT}/matrix$OCCUPANCY/${PREFIX}.partition | wc -l) loci." | tee -a ${DIR_OUTPUT}/summary.matrices

if [ ${ALIGN_TYPE} == "2" -a ${CODON12} == "1" ]; then
  echo "${OCCUPANCY}% occupancy matrix_cds12 has $(cat ${DIR_OUTPUT}/matrix${OCCUPANCY}/${PREFIX}.cds12.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat ${DIR_OUTPUT}/matrix$OCCUPANCY/${PREFIX}.cds12.partition | wc -l) loci." | tee -a ${DIR_OUTPUT}/summary.matrices

else
 echo
fi
