#!/bin/bash
#2021.08.04 wrote by ZF
#2023.03.29 revised by DSY

#Type 'bash BUSCO_extraction.sh'
#Tools TransDecoder and parallel are used in this script and will be automatically checked prior to formal analyses
#All the BUSCO results (run_* folders generated from BUSCO v3) are deposited in the same folder, e.g. BUSCOs/

##Checking the package dependency
echo "Checking the package dependency......"

#check TransDecoder
if [ $(which TransDecoder.LongOrfs) ]
    then
      echo "TransDecoder ...... OK"
      EXE_TRANSDECODER=$(which TransDecoder.LongOrfs)
      DIR_TRANSDECODER=${EXE_TRANSDECODER%/*}
    else
      until [ -x $DIR_TRANSDECODER/TransDecoder.LongOrfs ]
        do
          read -p "TransDecoder is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/TransDecoder-v5.5.0):      " DIR_TRANSDECODER_TEMP
          DIR_TRANSDECODER=$(realpath $(echo $DIR_TRANSDECODER_TEMP | sed "s/'//g"))
        done
      echo "TransDecoder ...... OK"
fi


#check parallel
if [ $(which parallel) ]
    then
      echo "parallel ...... OK"
      EXE_PARALLEL=$(which parallel)
      DIR_PARALLEL=${EXE_PARALLEL%/*}
    else
      until [ -x $DIR_PARALLEL/parallel ]
        do
          read -p "parallel is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/TransDecoder-v5.5.0):      " DIR_PARALLEL_TEMP
          DIR_PARALLEL=$(realpath $(echo $DIR_PARALLEL_TEMP | sed "s/'//g"))
        done
      echo "PARALLEL ...... OK"
fi


#Check the threads can be used
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
  until [ $THREADS -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done

#Input the version of BUSCO for assessments
read -p "Please input the version of BUSCO for assessments: 1. BUSCOv3; 2. BUSCOv5(version 5.4)   " BUSCO_VERSION
  until [ $BUSCO_VERSION -gt 0 -a $BUSCO_VERSION -lt 3 ]
    do
      read -p "Please input the version of BUSCO for assessments: 1. BUSCOv3; 2. BUSCOv5(version 5.4)      " BUSCO_VERSION
    done

#input the name of input directory
read -p "Please input the name of input directory containing all BUSCO results folders, e.g. /home/zf/Desktop/materials/datasets_examples/BUSCOs:      " DIR_INPUT_TEMP
DIR_INPUT=$(realpath $(echo $DIR_INPUT_TEMP | sed "s/'//g"))


if [ "$BUSCO_VERSION" == "1" ]; then
#copy the busco results to 0-raw_busco/
cp -r $DIR_INPUT 0-raw_busco

#Generate a file containing species list
ls 0-raw_busco/ | sed 's/run_//g' > species.list

#Generate a file of loci list
SPECIES_NAME=$(cat species.list)
for SPECIES in $SPECIES_NAME
do
cd 0-raw_busco/run_$SPECIES/single_copy_busco_sequences
ls | sed 's/.fna\|.faa//g' >> ../../../temp.list; cd ../../..
done
sort -n temp.list | uniq > loci.list
rm temp.list

#Modify the head name of the fasta files for each locus
export DIR_TRANSDECODER
TRANSDECODER_fun() {
	$DIR_TRANSDECODER/TransDecoder.LongOrfs -t $1 -m 20
        cp $1.transdecoder_dir/longest_orfs.cds $1
        cp $1.transdecoder_dir/longest_orfs.pep $(ls $1 | cut -d . -f1).faa
}
export -f TRANSDECODER_fun

for SPECIES in $SPECIES_NAME
	do
	  cd 0-raw_busco/run_$SPECIES/single_copy_busco_sequences/
	  sed -i -n '1,2p' *.f*a
          find . -name "*fna" | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 1 TRANSDECODER_fun %
          rm -rf *dir* *cmds
	  sed -i -n '1,2p' *.f*a
	  sed -i "1c >$SPECIES" *.f*a
	  sed -i "s/*//g" *.faa
	  cd ../../..
	done

#Merge sequences of the same locus into the fasta files
mkdir -p 1-raw_loci/fna 1-raw_loci/faa
LOCI_NAME=$(cat loci.list)
for LOCI in $LOCI_NAME; do touch 1-raw_loci/fna/$LOCI.fna && touch 1-raw_loci/faa/$LOCI.faa; done
for LOCI in $LOCI_NAME
do
    for SPECIES in $SPECIES_NAME
    do
    if [ -f 0-raw_busco/run_$SPECIES/single_copy_busco_sequences/$LOCI.fna ]
          then
            cat 0-raw_busco/run_$SPECIES/single_copy_busco_sequences/$LOCI.fna >> 1-raw_loci/fna/$LOCI.fna 
            cat 0-raw_busco/run_$SPECIES/single_copy_busco_sequences/$LOCI.faa >> 1-raw_loci/faa/$LOCI.faa
          else
            echo "$LOCI in $SPECIES does not exist" | tee -a 1-raw_loci/log.txt
       fi
    done
done

elif [ "$BUSCO_VERSION" == "2" ]; then
#copy the busco results to 0-raw_busco/
cp -r $DIR_INPUT 0-raw_busco

#Generate a file containing species list
ls 0-raw_busco/ | sed 's/ //g' > species.list

#Generate a file of loci list
SPECIES_NAME=$(cat species.list)
for SPECIES in $SPECIES_NAME
do
cd 0-raw_busco/$SPECIES/run_*/busco_sequences/single_copy_busco_sequences/
ls | sed 's/.fna\|.faa\|.gff//g' >> ../../../../../temp.list; cd ../../../../..
done
sort -n temp.list | uniq > loci.list
rm temp.list

#Modify the head name of the fasta files for each locus
export DIR_TRANSDECODER
TRANSDECODER_fun() {
  $DIR_TRANSDECODER/TransDecoder.LongOrfs -t $1 -m 20
        cp $1.transdecoder_dir/longest_orfs.cds $1
        cp $1.transdecoder_dir/longest_orfs.pep $(ls $1 | cut -d . -f1).faa
}
export -f TRANSDECODER_fun

for SPECIES in $SPECIES_NAME
  do
    cd 0-raw_busco/$SPECIES/run_*/busco_sequences/single_copy_busco_sequences/
    sed -i -n '1,2p' *.f*a
          find . -name "*fna" | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 1 TRANSDECODER_fun %
          rm -rf *dir* *cmds
    sed -i -n '1,2p' *.f*a
    sed -i "1c >$SPECIES" *.f*a
    sed -i "s/*//g" *.faa
    cd ../../../../..
  done

#Merge sequences of the same locus into the fasta files
mkdir -p 1-raw_loci/fna 1-raw_loci/faa
LOCI_NAME=$(cat loci.list)
for LOCI in $LOCI_NAME; do touch 1-raw_loci/fna/$LOCI.fna && touch 1-raw_loci/faa/$LOCI.faa; done
for LOCI in $LOCI_NAME
do
    for SPECIES in $SPECIES_NAME
    do
	  if [ -f 0-raw_busco/$SPECIES/run_*/busco_sequences/single_copy_busco_sequences/$LOCI.fna ]
          then
            cat 0-raw_busco/$SPECIES/run_*/busco_sequences/single_copy_busco_sequences/$LOCI.fna >> 1-raw_loci/fna/$LOCI.fna 
            cat 0-raw_busco/$SPECIES/run_*/busco_sequences/single_copy_busco_sequences/$LOCI.faa >> 1-raw_loci/faa/$LOCI.faa
          else
            echo "$LOCI in $SPECIES does not exist" | tee -a 1-raw_loci/log.txt
       fi
    done
done

fi

#Filter loci having too few taxa (less than three)
mkdir -p 2-loci_filter/fna 2-loci_filter/faa
TOTAL_TAXA=$(cat species.list | wc -l)
for LOCI in $LOCI_NAME
do
    echo -e $LOCI'\t'$(grep -o "^$LOCI" 1-raw_loci/log.txt | wc -l) | tee -a 2-loci_filter/sequence_number.log
    if [ $(grep -o "^$LOCI" 1-raw_loci/log.txt | wc -l) -lt `expr $TOTAL_TAXA - 2` ]
      then
        cp 1-raw_loci/fna/$LOCI.fna 2-loci_filter/fna
        cp 1-raw_loci/faa/$LOCI.faa 2-loci_filter/faa
        echo "$LOCI" >> 2-loci_filter/loci_name_filter.log
    fi
done
