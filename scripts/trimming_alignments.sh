#!/bin/bash
#2021.08.04 by ZF

#Trim alignments using trimal, BMGE and ClipKIT
#Type "bash trmming_alignments.sh"
#parallel and at least one of the three trimming tools should be installed
#ClipKIT were install by contructing CONDA virtual environment 'clipkit' and activated by 'source activate clipkit'


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


#Check the type of input alignments
read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide   " ALIGN_TYPE
  until [ $ALIGN_TYPE -gt 0 -a $ALIGN_TYPE -lt 3 ]
    do
      read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide      " ALIGN_TYPE
    done


#Check the trimming tool can be used
read -p "Please input the option for trimming tool: 1. trimAl; 2. BMGE; 3. ClipKIT   " TRIM_METHOD
  until [ $TRIM_METHOD -gt 0 -a $TRIM_METHOD -lt 4 ]
    do
      read -p "Please input the option for trimming tool: 1. trimAl; 2. BMGE; 3. ClipKIT      " TRIM_METHOD
    done


#Check the trimal installation and the method used by trimal
if [ "$TRIM_METHOD" == "1" ]; then
  #check trimAL installation
  if [ $(which trimal) ]
    then
      echo "trimal ...... OK"
      EXE_TRIMAL=$(which trimal)
      DIR_TRIMAL=${EXE_TRIMAL%/*}
    else
      until [ -x $DIR_TRIMAL/trimal ]
        do
          read -p "trimal is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/trimal-1.4.1/source):      " DIR_TRIMAL_TEMP
          DIR_TRIMAL=$(realpath $(echo $DIR_TRIMAL_TEMP | sed "s/'//g"))
        done
      echo "trimal ...... OK"
  fi

  read -p "Please input the option for trimming strategy: 1. -automated1; 2. -gappyout   " TRIMAL_METHOD
    until [ $TRIMAL_METHOD -gt 0 -a $TRIMAL_METHOD -lt 3 ]
      do
        read -p "Please input the option for trimming strategy: 1. -automated1; 2. -gappyout      " TRIMAL_METHOD
      done
fi


#Check the BMGE installation and the method used by BMGE
if [ "$TRIM_METHOD" == "2" ]; then
  #check BMGE installation
  if [ $(which BMGE.jar) ]
    then
      echo "BMGE ...... OK"
      EXE_BMGE=$(which BMGE.jar)
      DIR_BMGE=${EXE_BMGE%/*}
    else
      until [ -s $DIR_BMGE/BMGE.jar ]
        do
          read -p "BMGE is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/BMGE-1.12):      " DIR_BMGE_TEMP
          DIR_BMGE=$(realpath $(echo $DIR_BMGE_TEMP | sed "s/'//g"))
        done
      echo "BMGE ...... OK"
  fi

  read -p "Please input the option for trimming strategy: 1. default; 2. stringent (-h 0.4)   " BMGE_METHOD
    until [ $BMGE_METHOD -gt 0 -a $BMGE_METHOD -lt 3 ]
      do
        read -p "Please input the option for trimming strategy: 1. default; 2. stringent (-h 0.4)      " BMGE_METHOD
      done
fi


#input the name of output directory
read -p "Please input the name of output directory:      " DIR_OUTPUT_TEMP
DIR_OUTPUT=$(echo $DIR_OUTPUT_TEMP | sed "s/'//g")
test -d $DIR_OUTPUT || mkdir $DIR_OUTPUT


#Check the ClipKIT installation and the method used by ClipKIT
if [ "$TRIM_METHOD" == "3" ]; then
  #check ClipKIT installation
  if [ $(which clipkit) ]
    then
      echo "clipkit ...... OK"
      EXE_CLIPKIT=$(which clipkit)
      DIR_CLIPKIT=${EXE_CLIPKIT%/*}
    else
      until [ -x $DIR_CLIPKIT/clipkit ]
        do
          read -p "ClipKIT is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/ClipKIT/):      " DIR_CLIPKIT_TEMP
          DIR_CLIPKIT=$(realpath $(echo $DIR_CLIPKIT_TEMP | sed "s/'//g"))
        done
      echo "ClipKIT ...... OK"
  fi

  #check the method used by ClipKIT
  read -p "Please input the option for trimming strategy: 1. smart-gap (dynamic determination of gaps threshold); 2. kpi-smart-gap (a combination of keeping parsimony informative sites and smart-gap-based trimming); 3. kpic-smart-gap (a combination of keeping parismony informative and constant sites and smart-gap-based trimming   " CLIPKIT_METHOD
    until [ $CLIPKIT_METHOD -gt 0 -a $CLIPKIT_METHOD -lt 4 ]
      do
        read -p "Please input the option for trimming strategy: 1. smart-gap (dynamic determination of gaps threshold); 2. kpi-smart-gap (a combination of keeping parsimony informative sites and smart-gap-based trimming); 3. kpic-smart-gap (a combination of keeping parismony informative and constant sites and smart-gap-based trimming      " CLIPKIT_METHOD
      done

  #ClipKIT trimming
  read -p "Please input the folder names (with its path) containing aminoacid/nucleotide alignments, for example '3-faa_align':      " INPUT_TEMP
  INPUT=$(realpath $(echo $INPUT_TEMP | sed "s/'//g"))
  [ "$ALIGN_TYPE" == "1" ] && mkdir -p $DIR_OUTPUT/clipkit-faa || mkdir -p $DIR_OUTPUT/clipkit-fna
  ls $INPUT/ | cut -d "." -f1 > $DIR_OUTPUT/loci.list
  for loci in $(cat $DIR_OUTPUT/loci.list); do cp $INPUT/$loci.* $DIR_OUTPUT/$loci; done
  cd $DIR_OUTPUT
  if [ "$CLIPKIT_METHOD" == "1" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_CLIPKIT/clipkit {} -o {}.fas -m smart-gap
  elif [ "$CLIPKIT_METHOD" == "2" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_CLIPKIT/clipkit {} -o {}.fas -m kpi-smart-gap
  elif [ "$CLIPKIT_METHOD" == "3" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_CLIPKIT/clipkit {} -o {}.fas -m kpic-smart-gap
  fi
  cat loci.list | $DIR_PARALLEL/parallel -j $THREADS rm {}
  [ "$ALIGN_TYPE" == "1" ] && mv *fas clipkit-faa/ || mv *fas clipkit-fna/

for loci in $(cat loci.list)
do 
    a=$(cat clipkit-faa/$loci.fas | sed -n 2p) 
    b=$(echo ${a:0:1})
    if [ "$b" == ">" ] 
    then
     echo $loci >> missing.loci
    fi
done

if [ -f missing.loci ]
then
mkdir missing_data
for missing in $(cat missing.loci)
do
   mv clipkit-faa/$missing.fas missing_data
done
ls clipkit-faa/ |  sed "s/.fas//g" > loci_reserve.list
fi

fi




if [ "$TRIM_METHOD" == "1" ] && [ "$ALIGN_TYPE" == "1" ]; then
  read -p "Please input the folder names (with its path) containing amino acid alignments and unaligned nucleotide sequences, for example '3-faa_align 2-loci_filter/fna':      " INPUT_TEMP1 INPUT_TEMP2
  INPUT_AA=$(realpath $(echo $INPUT_TEMP1 | sed "s/'//g"))
  INPUT_NUC=$(realpath $(echo $INPUT_TEMP2 | sed "s/'//g"))
  mkdir -p $DIR_OUTPUT/trimal_faa $DIR_OUTPUT/trimal_fna
  cp $INPUT_AA/* $DIR_OUTPUT
  ls $INPUT_AA/ | cut -d "." -f1 > $DIR_OUTPUT/loci.list
  for loci in $(cat $DIR_OUTPUT/loci.list); do cp $INPUT_AA/$loci.* $DIR_OUTPUT/$loci.faa; done
  for loci in $(cat $DIR_OUTPUT/loci.list); do cp $INPUT_NUC/$loci.* $DIR_OUTPUT/$loci.fna; done
  cd $DIR_OUTPUT
  if [ "$TRIMAL_METHOD" == "1" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_TRIMAL/trimal -in {}.faa -out {}.fas -automated1
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_TRIMAL/trimal -in {}.faa -out {}.nuc.fas -automated1 -backtrans {}.fna
  elif [ "$TRIMAL_METHOD" == "2" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_TRIMAL/trimal -in {}.faa -out {}.fas -gappyout
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_TRIMAL/trimal -in {}.faa -out {}.nuc.fas -gappyout -backtrans {}.fna
  fi
  cat loci.list | $DIR_PARALLEL/parallel -j $THREADS rm {}.faa {}.fna
  mv *nuc.fas trimal_fna && mv *fas trimal_faa

for loci in $(cat loci.list)
do 
    a=$(cat trimal_faa/$loci.fas | sed -n 2p) 
    b=$(echo ${a:0:1})
    if [ "$b" == ">" ] 
    then
     echo $loci >> missing.loci
    fi
done

if [ -f missing.loci ]
then
mkdir missing-faa_data missing-fna_data
for missing in $(cat missing.loci)
do
   mv trimal_faa/$missing.fas missing-faa_data
   mv trimal_fna/$missing.nuc.fas missing-fna_data
done
ls trimal_faa/ |  sed "s/.fas//g" > loci_reserve.list
fi

fi




if [ "$TRIM_METHOD" == "1" ] && [ "$ALIGN_TYPE" == "2" ]; then
  read -p "Please input the folder names (with its path) containing nucleotide alignments, for example '2-loci_filter/fna':      " INPUT_TEMP
  INPUT_NUC=$(realpath $(echo $INPUT_TEMP | sed "s/'//g"))
  mkdir -p $DIR_OUTPUT/trimal_fna
  cp $INPUT_NUC/* $DIR_OUTPUT
  ls $INPUT_NUC/ | cut -d "." -f1 > $DIR_OUTPUT/loci.list
  for loci in $(cat $DIR_OUTPUT/loci.list); do cp $INPUT_NUC/$loci.* $DIR_OUTPUT/$loci.fna; done
  cd $DIR_OUTPUT
  if [ "$TRIMAL_METHOD" == "1" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_TRIMAL/trimal -in {}.fna -out {}.fas -automated1
  elif [ "$TRIMAL_METHOD" == "2" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_TRIMAL/trimal -in {}.fna -out {}.fas -gappyout
  fi
  cat loci.list | $DIR_PARALLEL/parallel -j $THREADS rm {}.fna
  mv *fas trimal_fna

for loci in $(cat loci.list)
do 
    a=$(cat trimal_fna/$loci.nuc.fas | sed -n 2p) 
    b=$(echo ${a:0:1})
    if [ "$b" == ">" ] 
    then
     echo $loci >> missing.loci
    fi
done

if [ -f missing.loci ]
then
mkdir missing_data
for missing in $(cat missing.loci)
do
   mv trimal_fna/$missing.nuc.fas missing_data
done
ls trimal_fna/ |  sed "s/.nuc.fas//g" > loci_reserve.list
fi

fi


if [ "$TRIM_METHOD" == "2" ] && [ "$ALIGN_TYPE" == "1" ]; then
  read -p "Please input the folder names (with its path) containing amino acid alignments, for example '3-faa_align':      " INPUT_TEMP
  INPUT=$(realpath $(echo $INPUT_TEMP | sed "s/'//g"))
  mkdir -p $DIR_OUTPUT/bmge_faa
  ls $INPUT/ | cut -d "." -f1 > $DIR_OUTPUT/loci.list
  for loci in $(cat $DIR_OUTPUT/loci.list); do cp $INPUT/$loci.* $DIR_OUTPUT/$loci; done
  cd $DIR_OUTPUT
  if [ "$BMGE_METHOD" == "1" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS java -jar $DIR_BMGE/BMGE.jar -i {} -t AA -of {}.fas
  elif [ "$BMGE_METHOD" == "2" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS java -jar $DIR_BMGE/BMGE.jar -i {} -t AA -m BLOSUM90 -h 0.4 -of {}.fas
  fi
  cat loci.list | $DIR_PARALLEL/parallel -j $THREADS rm {}
  mv *fas bmge_faa

for loci in $(cat loci.list)
do 
    a=$(cat bmge_faa/$loci.fas | sed -n 2p) 
    b=$(echo ${a:0:1})
    if [ "$b" == ">" ] 
    then
     echo $loci >> missing.loci
    fi
done

if [ -f missing.loci ]
then
mkdir missing_data
for missing in $(cat missing.loci)
do
   mv bmge_faa/$missing.fas missing_data
done
ls bmge_faa/ |  sed "s/.fas//g" > loci_reserve.list
fi

fi


if [ "$TRIM_METHOD" == "2" ] && [ "$ALIGN_TYPE" == "2" ]; then
  read -p "Please input the folder names (with its path) containing nucleotide alignments, for example '3-faa_align':      " INPUT_TEMP
  INPUT=$(realpath $(echo $INPUT_TEMP | sed "s/'//g"))
  mkdir -p $DIR_OUTPUT/bmge_fna
  ls $INPUT/ | cut -d "." -f1 > $DIR_OUTPUT/loci.list
  for loci in $(cat $DIR_OUTPUT/loci.list); do cp $INPUT/$loci.* $DIR_OUTPUT/$loci; done
  cd $DIR_OUTPUT
  if [ "$BMGE_METHOD" == "1" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS java -jar $DIR_BMGE/BMGE.jar -i {} -t DNA -of {}.fas
  elif [ "$BMGE_METHOD" == "2" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS java -jar $DIR_BMGE/BMGE.jar -i {} -t DNA -m DNAPAM50 -h 0.4 -of {}.fas
  fi
  cat loci.list | $DIR_PARALLEL/parallel -j $THREADS rm {}
  mv *fas bmge_fna

for loci in $(cat loci.list)
do 
    a=$(cat bmge_fna/$loci.fas | sed -n 2p) 
    b=$(echo ${a:0:1})
    if [ "$b" == ">" ] 
    then
     echo $loci >> missing.loci
    fi
done

if [ -f missing.loci ]
then
mkdir missing_data
for missing in $(cat missing.loci)
do
   mv bmge_fna/$missing.fas missing_data
done
ls bmge_fna/ |  sed "s/.fas//g" > loci_reserve.list
fi

fi
