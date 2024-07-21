#!/bin/bash
#2021.08.04 wrote by ZF
#2023.06.29 revised by DSY

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
  read -p "Please input the folder names (with its path) containing aminoacid/nucleotide alignments, for example '3-faa_align/faa':      " DIR_INPUT_TEMP
  DIR_INPUT=$(realpath $(echo $DIR_INPUT_TEMP | sed "s/'//g"))
  if [ "$ALIGN_TYPE" == "1" ]; then
  mkdir -p $DIR_OUTPUT/clipkit-faa 
  mkdir -p $DIR_OUTPUT/clipkit-fna
  elif [ "$ALIGN_TYPE" == "2" ]; then
  mkdir -p $DIR_OUTPUT/clipkit-fna
  fi
  ls $DIR_INPUT/ > $DIR_OUTPUT/loci.list
  cp -r $DIR_INPUT $DIR_OUTPUT/temp_clipkit && cd $DIR_OUTPUT/temp_clipkit
  if [ "$CLIPKIT_METHOD" == "1" ]; then
    cat ../loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_CLIPKIT/clipkit {} -o ../{} -l -m smart-gap
  elif [ "$CLIPKIT_METHOD" == "2" ]; then
    cat ../loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_CLIPKIT/clipkit {} -o ../{} -l -m kpi-smart-gap
  elif [ "$CLIPKIT_METHOD" == "3" ]; then
    cat ../loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_CLIPKIT/clipkit {} -o ../{} -l -m kpic-smart-gap
  fi
  cd ../
  rm -rf temp_clipkit/

#Trim the the corresponding nucleotide sequences
if [ "$ALIGN_TYPE" == "1" ]; then
echo
read -p "Do you want to trim the the corresponding nucleotide sequences? 1. Yes; 2. No:      " NUC_TRIM
  until [ $NUC_TRIM -gt 0 -a $NUC_TRIM -lt 3 ]
    do
      read -p "Do you want to trim the the corresponding nucleotide sequences? 1. Yes; 2. No:      " NUC_TRIM
    done

  if [ "$NUC_TRIM" == "1" ]; then
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
   read -p "Please input the folder names (with its path) containing nucleotide alignments, for example '3-faa_align/fna':      " INPUT_TEMP
   INPUT_NUC=$(realpath $(echo $INPUT_TEMP | sed "s/'//g"))
   for LOCI in $(cat loci.list); do mv $LOCI clipkit-faa/; done
   mkdir temp_nuc
   cp -r $INPUT_NUC/* temp_nuc
   TRIM_fun() {
    pk_thread_dna -p clipkit-faa/$1 -n temp_nuc/$1 -c $1.log > clipkit-fna/$1
  }
  export -f TRIM_fun
  cat loci.list | $DIR_PARALLEL/parallel -j $THREADS TRIM_fun
  rm -rf temp_nuc/ *.log
  elif [ "$NUC_TRIM" == "2" ]; then
    for LOCI in $(cat loci.list); do mv $LOCI clipkit-faa/; done
  rm -rf *.log
  fi
  elif [ "$ALIGN_TYPE" == "2" ]; then
  for LOCI in $(cat loci.list); do mv $LOCI clipkit-fna/; done
  rm -rf *.log
fi

if [ "$ALIGN_TYPE" == "1" ]; then
  
  if [ "$NUC_TRIM" == "1" ]; then
   for loci in $(cat loci.list)
   do 
     a=$(cat clipkit-faa/$loci | sed -n 2p) 
     b=$(echo ${a:0:1})
     if [ "$b" == ">" ] 
     then
       echo $loci >> missing.loci
     fi
   done
   if [ -f missing.loci ]
     then
     for missing in $(cat missing.loci)
     do
     rm clipkit-faa/$missing 
     rm clipkit-fna/$missing 
     done
     rm loci.list
     ls clipkit-faa/ > loci.list
   fi
  elif [ "$NUC_TRIM" == "2" ]; then
    for loci in $(cat loci.list)
   do 
     a=$(cat clipkit-faa/$loci | sed -n 2p) 
     b=$(echo ${a:0:1})
     if [ "$b" == ">" ] 
     then
       echo $loci >> missing.loci
     fi
   done
   if [ -f missing.loci ]
     then
     for missing in $(cat missing.loci)
     do
     rm clipkit-faa/$missing 
     done
     rm -rf loci.list clipkit-fna/
     ls clipkit-faa/ > loci.list
   fi
  fi
elif [ "$ALIGN_TYPE" == "2" ]; then
for loci in $(cat loci.list)
do 
    a=$(cat clipkit-fna/$loci | sed -n 2p) 
    b=$(echo ${a:0:1})
    if [ "$b" == ">" ] 
    then
     echo $loci >> missing.loci
    fi
done
if [ -f missing.loci ]
then
for missing in $(cat missing.loci)
do
   rm clipkit-fna/$missing
done
rm loci.list
ls clipkit-faa/ > loci.list
fi
fi
fi

if [ "$TRIM_METHOD" == "1" ] && [ "$ALIGN_TYPE" == "1" ]; then
  read -p "Please input the folder names (with its path) containing amino acid alignments and unaligned nucleotide sequences, for example '3-faa_align 2-loci_filter/fna':      " INPUT_TEMP1 INPUT_TEMP2
  INPUT_AA=$(realpath $(echo $INPUT_TEMP1 | sed "s/'//g"))
  INPUT_NUC=$(realpath $(echo $INPUT_TEMP2 | sed "s/'//g"))
  mkdir -p $DIR_OUTPUT/trimal_faa $DIR_OUTPUT/trimal_fna $DIR_OUTPUT/temp1 $DIR_OUTPUT/temp2
  cd $DIR_OUTPUT
  cp $INPUT_AA/* temp1/
  cp $INPUT_NUC/* temp2/
  rename "s/.fna/.fas/" temp2/*
  ls temp1/ > loci.list
  if [ "$TRIMAL_METHOD" == "1" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_TRIMAL/trimal -in temp1/{} -out trimal_faa/{} -automated1
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_TRIMAL/trimal -in temp1/{} -out trimal_fna/{} -automated1 -backtrans temp2/{}
  elif [ "$TRIMAL_METHOD" == "2" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_TRIMAL/trimal -in temp1/{} -out trimal_faa/{} -gappyout
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_TRIMAL/trimal -in temp1/{} -out trimal_fna/{} -gappyout -backtrans temp2/{}
  fi
  rm -rf temp1/ temp2/

for loci in $(cat loci.list)
do 
    a=$(cat trimal_faa/$loci | sed -n 2p) 
    b=$(echo ${a:0:1})
    if [ "$b" == ">" ] 
    then
     echo $loci >> missing.loci
    fi
done

if [ -f missing.loci ]
then
for missing in $(cat missing.loci)
do
   rm trimal_faa/$missing
   rm trimal_fna/$missing
done
rm loci.list
ls trimal_faa/ > loci.list
fi
fi


if [ "$TRIM_METHOD" == "1" ] && [ "$ALIGN_TYPE" == "2" ]; then
  read -p "Please input the folder names (with its path) containing nucleotide alignments, for example '2-loci_filter/fna':      " INPUT_TEMP
  INPUT_NUC=$(realpath $(echo $INPUT_TEMP | sed "s/'//g"))
  mkdir -p $DIR_OUTPUT/trimal_fna $DIR_OUTPUT/temp
  cd $DIR_OUTPUT
  cp $INPUT_NUC/* temp/
  ls temp/ > loci.list
  if [ "$TRIMAL_METHOD" == "1" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_TRIMAL/trimal -in temp/{} -out trimal_fna/{} -automated1
  elif [ "$TRIMAL_METHOD" == "2" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_TRIMAL/trimal -in temp/{} -out trimal_fna/{} -gappyout
  fi
  rm -rf temp/

for loci in $(cat loci.list)
do 
    a=$(cat trimal_fna/$loci | sed -n 2p) 
    b=$(echo ${a:0:1})
    if [ "$b" == ">" ] 
    then
     echo $loci >> missing.loci
    fi
done

if [ -f missing.loci ]
then
for missing in $(cat missing.loci)
do
   rm trimal_fna/$missing
done
rm loci.list
ls trimal_fna/ > loci.list
fi
fi


if [ "$TRIM_METHOD" == "2" ] && [ "$ALIGN_TYPE" == "1" ]; then
  read -p "Please input the folder names (with its path) containing amino acid alignments, for example '3-faa_align':      " INPUT_TEMP
  INPUT=$(realpath $(echo $INPUT_TEMP | sed "s/'//g"))
  mkdir -p $DIR_OUTPUT/bmge_faa
  ls $INPUT/ > $DIR_OUTPUT/loci.list
  for loci in $(cat $DIR_OUTPUT/loci.list); do cp $INPUT/$loci $DIR_OUTPUT/; done
  cd $DIR_OUTPUT
  if [ "$BMGE_METHOD" == "1" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS java -jar $DIR_BMGE/BMGE.jar -i {} -t AA -of bmge_faa/{}
  elif [ "$BMGE_METHOD" == "2" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS java -jar $DIR_BMGE/BMGE.jar -i {} -t AA -m BLOSUM90 -h 0.4 -of bmge_faa/{}
  fi
  cat loci.list | $DIR_PARALLEL/parallel -j $THREADS rm {}

for loci in $(cat loci.list)
do 
    a=$(cat bmge_faa/$loci | sed -n 2p) 
    b=$(echo ${a:0:1})
    if [ "$b" == ">" ] 
    then
     echo $loci >> missing.loci
    fi
done

if [ -f missing.loci ]
then
for missing in $(cat missing.loci)
do
   rm bmge_faa/$missing
done
rm loci.list
ls bmge_faa/ > loci.list
fi
fi


if [ "$TRIM_METHOD" == "2" ] && [ "$ALIGN_TYPE" == "2" ]; then
  read -p "Please input the folder names (with its path) containing nucleotide alignments, for example '3-faa_align':      " INPUT_TEMP
  INPUT=$(realpath $(echo $INPUT_TEMP | sed "s/'//g"))
  mkdir -p $DIR_OUTPUT/bmge_fna
  ls $INPUT/ > $DIR_OUTPUT/loci.list
  for loci in $(cat $DIR_OUTPUT/loci.list); do cp $INPUT/$loci $DIR_OUTPUT/; done
  cd $DIR_OUTPUT
  if [ "$BMGE_METHOD" == "1" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS java -jar $DIR_BMGE/BMGE.jar -i {} -t DNA -of bmge_fna/{}
  elif [ "$BMGE_METHOD" == "2" ]; then
    cat loci.list | $DIR_PARALLEL/parallel -j $THREADS java -jar $DIR_BMGE/BMGE.jar -i {} -t DNA -m DNAPAM50 -h 0.4 -of bmge_fna/{}
  fi
  cat loci.list | $DIR_PARALLEL/parallel -j $THREADS rm {}

for loci in $(cat loci.list)
do 
    a=$(cat bmge_fna/$loci | sed -n 2p) 
    b=$(echo ${a:0:1})
    if [ "$b" == ">" ] 
    then
     echo $loci >> missing.loci
    fi
done

if [ -f missing.loci ]
then
for missing in $(cat missing.loci)
do
   rm bmge_fna/$missing
done
rm loci.list
ls bmge_fna/ > loci.list
fi
fi
