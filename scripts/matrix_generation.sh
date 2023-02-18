#!/bin/bash
#2021.08.04 wrote by ZF, 2023.1.11 modified by DSY

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
      until [ -x $DIR_PHYKIT/phykit ]
        do
          read -p "PhyKIT is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_PHYKIT_TEMP
          DIR_PHYKIT=$(realpath $(echo $DIR_PHYKIT_TEMP | sed "s/'//g"))
        done
      echo "PhyKIT ...... OK"
fi


DIR_CURR=$(echo $PWD)

#input the name of input directory
read -p "Please input the name of input directory containing all alignments, e.g. 4-trim/clipkit-kpi:      " DIR_INPUT_TEMP
DIR_INPUT=$(realpath $(echo $DIR_INPUT_TEMP | sed "s/'//g"))


#input the name of output directory
read -p "Please input the name of output directory, or an existing directory:      " DIR_OUTPUT_TEMP
DIR_OUTPUT_TEMP1=$(echo $DIR_OUTPUT_TEMP | sed "s/'//g")
test -d $DIR_OUTPUT_TEMP1 || mkdir -p $DIR_OUTPUT_TEMP1
cd $DIR_OUTPUT_TEMP1 && DIR_OUTPUT=$(echo $PWD) && cd $DIR_CURR


#input the prefix for the generated matrix-related files
read -p "Please input the name of a prefix for the generated matrix-related files, e.g. DATASET1:      " PREFIX


#input the taxon occupancy values for the matrices
read -p "Please input the minimum percentage value for taxa occupancy, usually ranging from 50% to 100%, e.g. 50, 75, 90:      " OCCUPANCY


#Check the type of input alignments
read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide   " ALIGN_TYPE
  until [ $ALIGN_TYPE -gt 0 -a $ALIGN_TYPE -lt 3 ]
    do
      read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide      " ALIGN_TYPE
    done


if [ "$ALIGN_TYPE" == "1" ]; then
cd $DIR_INPUT
ls > ../$PREFIX.alignments
$DIR_PHYKIT/phykit cc -a ../$PREFIX.alignments -p $PREFIX
rm ../$PREFIX.alignments

mkdir -p $DIR_OUTPUT/all $DIR_OUTPUT/matrix$OCCUPANCY/alignments
mv $PREFIX.* $DIR_OUTPUT/all

cd $DIR_OUTPUT/all 
TOTAL_LINE=$(cat $PREFIX.occupancy | wc -l)
for line in $(seq $TOTAL_LINE)
  do
    loci=$(sed -n "$line"p $PREFIX.occupancy | awk '{print $1}')
    num=$(sed -n "$line"p $PREFIX.occupancy | awk '{print $4}')
    OCCU=$(echo "scale=4;($OCCUPANCY/100)"|bc)
    diff=$(echo "scale=4;($num-$OCCU)"|bc)
    num1=`echo "$diff < 0" |bc`
    test "$num1" = 0 && echo $loci >> $DIR_OUTPUT/matrix$OCCUPANCY/$PREFIX.alignments
  done

cd $DIR_OUTPUT/matrix$OCCUPANCY/alignments
for file in $(cat ../$PREFIX.alignments); do cp $DIR_INPUT/$file ./; done
$DIR_PHYKIT/phykit cc -a ../$PREFIX.alignments -p $PREFIX
mv  $PREFIX.* $DIR_OUTPUT/matrix$OCCUPANCY/
rm -rf $DIR_OUTPUT/all

elif [ "$ALIGN_TYPE" == "2" ]; then
#check FASconCAT-G_v1.04.pl
if [ $(which FASconCAT-G_v1.04.pl) ]
    then
      echo "FASconCAT-G_v1.04.pl ...... OK"
      EXE_FASCONCAT=$(which FASconCAT-G_v1.04.pl)
      DIR_FASCONCAT=${EXE_FASCONCAT%/*}
    else
      until [ -x $DIR_FASCONCAT/FASconCAT-G_v1.04.pl ]
        do
          read -p "FASconCAT-G_v1.04.pl is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/FASconCAT-G-1.04/FASconCAT-G_v1.04.pl):      " DIR_FASCONCAT_TEMP
          DIR_FASCONCAT=$(realpath $(echo $DIR_FASCONCAT_TEMP | sed "s/'//g"))
        done
      echo "FASconCAT-G_v1.04.pl ...... OK"
fi

  mkdir -p $DIR_OUTPUT/all/alignments
  cd $DIR_OUTPUT/all/alignments
  cp $DIR_INPUT/* ./
  ls > ../$PREFIX.alignments
  $DIR_PHYKIT/phykit cc -a ../$PREFIX.alignments -p $PREFIX
  rm ../$PREFIX.alignments
  mv $PREFIX* ..
  cd ..
  TOTAL_LINE=$(cat $PREFIX.occupancy | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p $PREFIX.occupancy | awk '{print $1}')
      num=$(sed -n "$line"p $PREFIX.occupancy | awk '{print $4}')
      OCCU=$(echo "scale=4;($OCCUPANCY/100)"|bc)
      diff=$(echo "scale=4;($num-$OCCU)"|bc)
      num1=`echo "$diff < 0" |bc`
      test "$num1" = 0 && echo $loci >> $DIR_OUTPUT/all/$PREFIX.alignments
    done

  mkdir -p $DIR_OUTPUT/cds12_fna/alignments $DIR_OUTPUT/cds_fna/alignments $DIR_OUTPUT/fas
  cd $DIR_OUTPUT/fas
  cp $DIR_INPUT/* .
  cp $DIR_OUTPUT/all/$PREFIX.alignments $DIR_OUTPUT/cds12_fna
  cp $DIR_OUTPUT/all/$PREFIX.alignments $DIR_OUTPUT/cds_fna

  cd $DIR_OUTPUT/cds_fna/alignments
  for LOCI in $(cat ../*.alignments)
  do
    cp ../../fas/$LOCI .
  done
  $DIR_PHYKIT/phykit cc -a ../*.alignments -p $PREFIX
  mv $PREFIX* ..
  rm ../*.occupancy
  cd $DIR_OUTPUT/

  cd $DIR_OUTPUT/cds12_fna/alignments
  for LOCI in $(cat ../*.alignments)
  do 
    cp ../../fas/$LOCI .
  done
  perl $DIR_FASCONCAT/FASconCAT-G*.pl -s -l -d
  mv *.txt ../$PREFIX.partition
  mv FcC_supermatrix.fas ../$PREFIX.fa
  rm *.xls
  cd $DIR_OUTPUT
  rm -rf fas/ all/
fi

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
      until [ -x $DIR_SEQKIT/seqkit ]
        do
          read -p "SeqKit is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_SEQKIT_TEMP
          DIR_SEQKIT=$(echo $DIR_SEQKIT_TEMP | sed "s/'//g")
        done
      echo "SeqKit ...... OK"
  fi

  read -p "Please input the species name of outgroup taxon, such as 'Zootermopsis_nevadensis':      "  OUTGROUP
  if [ "$ALIGN_TYPE" == "1" ]; then  
  cd $DIR_OUTPUT/matrix$OCCUPANCY/
  echo "$OUTGROUP" > list.outgroup
  cat $PREFIX.fa | $DIR_SEQKIT/seqkit grep -f list.outgroup | $DIR_SEQKIT/seqkit seq -u -w 0 > outgroup.fa
  cat $PREFIX.fa | $DIR_SEQKIT/seqkit grep -v -f list.outgroup | $DIR_SEQKIT/seqkit seq -u -w 0 > ingroup.fa
  cat outgroup.fa ingroup.fa > $PREFIX.fa
  rm ingroup.fa outgroup.fa list.outgroup
  echo
  elif [ "$ALIGN_TYPE" == "2" ]; then
  cd $DIR_OUTPUT/
  echo "$OUTGROUP" > list.outgroup
  cat cds_fna/$PREFIX.fa | $DIR_SEQKIT/seqkit grep -f list.outgroup | $DIR_SEQKIT/seqkit seq -u -w 0 > cds_fna/outgroup.fa
  cat cds_fna/$PREFIX.fa | $DIR_SEQKIT/seqkit grep -v -f list.outgroup | $DIR_SEQKIT/seqkit seq -u -w 0 > cds_fna/ingroup.fa
  cat cds_fna/outgroup.fa cds_fna/ingroup.fa > cds_fna/$PREFIX.fa
  cat cds12_fna/$PREFIX.fa | $DIR_SEQKIT/seqkit grep -f list.outgroup | $DIR_SEQKIT/seqkit seq -u -w 0 > cds12_fna/outgroup.fa
  cat cds12_fna/$PREFIX.fa | $DIR_SEQKIT/seqkit grep -v -f list.outgroup | $DIR_SEQKIT/seqkit seq -u -w 0 > cds12_fna/ingroup.fa
  cat cds12_fna/outgroup.fa cds12_fna/ingroup.fa > cds12_fna/$PREFIX.fa
  rm cds_fna/ingroup.fa cds_fna/outgroup.fa list.outgroup cds12_fna/ingroup.fa cds12_fna/outgroup.fa
  fi
else
  echo
fi

if [ "$ALIGN_TYPE" == "1" ]; then
echo "Individual loci alignments, concatenated matrix and partition file are deposited in the OUTPUT/matrix$OCCUPANCY."
echo "$OCCUPANCY% occupancy matrix has $(cat $DIR_OUTPUT/matrix$OCCUPANCY/$PREFIX.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/matrix$OCCUPANCY/$PREFIX.partition | wc -l) loci." | tee -a $DIR_OUTPUT/summary.matrices
cd $DIR_CURR
elif [ "$ALIGN_TYPE" == "2" ]; then
echo "Individual loci alignments, concatenated matrix and partition file are deposited in the OUTPUT/matrix$OCCUPANCY."
echo "$OCCUPANCY% occupancy matrix has $(cat $DIR_OUTPUT/cds_fna/$PREFIX.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/cds_fna/$PREFIX.partition | wc -l) loci with nucleotides." | tee -a $DIR_OUTPUT/summary.matrices
echo "$OCCUPANCY% occupancy matrix has $(cat $DIR_OUTPUT/cds12_fna/$PREFIX.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/cds12_fna/$PREFIX.partition | wc -l) loci with nucleotides reject 3rd codon position." | tee -a $DIR_OUTPUT/summary.matrices
cd $DIR_CURR
fi
