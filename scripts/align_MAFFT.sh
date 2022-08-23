#!/bin/bash
#2021.08.04 wrote by ZF
#2022.08.18 revised by DSY

#Type 'bash align_MAFFT.sh'
#Tools MAFFT and MAGUS may be used in this script and will be automatically checked prior to formal analyses

#check MAFFT
if [ $(which mafft) ]
    then
      echo "MAFFT ...... OK"
      EXE_MAFFT=$(which mafft)
      DIR_MAFFT=${EXE_MAFFT%/*}
    else
      until [ -x $DIR_MAFFT/mafft ]
        do
          read -p "MAFFT is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_MAFFT_TEMP
          DIR_MAFFT=$(realpath $(echo $DIR_MAFFT_TEMP | sed "s/'//g"))
        done
      echo "MAFFT ...... OK"
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


#Check the MAFFT method can be used
read -p "Please input the option for MAFFT-based strategy: 1. mafft-auto; 2. linsi; 3. einsi; 4. ginsi; 5. MAGUS   " MAFFT_METHOD
  until [ $MAFFT_METHOD -gt 0 -a $MAFFT_METHOD -lt 6 ]
    do
      read -p "Please input the option for MAFFT-based strategy: 1. mafft-auto; 2. linsi; 3. einsi; 4. ginsi; 5. MAGUS      " MAFFT_METHOD
    done


#Check the threads can be used
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
  until [ $THREADS -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done


#Check the type of input alignments
read -p "Please input the option for input unalignments: 1. amino acid; 2. nucleotide   " ALIGN_TYPE
  until [ $ALIGN_TYPE -gt 0 -a $ALIGN_TYPE -lt 3 ]
    do
      read -p "Please input the option for input unalignments: 1. amino acid; 2. nucleotide      " ALIGN_TYPE
    done


DIR_CURR=$(echo $PWD)

#input the name of input directory
read -p "Please input the name of input directory containing all unalignments, e.g. 2-loci_filter/faa:      " DIR_INPUT_TEMP
DIR_INPUT=$(realpath $(echo $DIR_INPUT_TEMP | sed "s/'//g"))


#input the name of output directory
read -p "Please input the name of output directory, or an existing directory:      " DIR_OUTPUT_TEMP
DIR_OUTPUT_TEMP1=$(echo $DIR_OUTPUT_TEMP | sed "s/'//g")
test -d $DIR_OUTPUT_TEMP1 || mkdir -p $DIR_OUTPUT_TEMP1
cd $DIR_OUTPUT_TEMP1 && DIR_OUTPUT=$(echo $PWD) && cd $DIR_CURR


TEMP=$(ls "$DIR_INPUT")


if [ "$MAFFT_METHOD" == "5" ]; then
  #check MAGUS
  if [ $(which magus.py) ]
    then
      echo "MAGUS ...... OK"
      EXE_MAGUS=$(which magus.py)
      DIR_MAGUS=${EXE_MAGUS%/*}
    else
      until [ -s $DIR_MAGUS/magus.py ]
        do
          read -p "MAGUS is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/MAGUS):      " DIR_MAGUS_TEMP
          DIR_MAGUS=$(realpath $(echo $DIR_MAGUS_TEMP | sed "s/'//g"))
        done
      echo "MAGUS ...... OK"
  fi

  cp -r $DIR_INPUT temp_magus && cd temp_magus
  for file in $TEMP
    do
      python3 $DIR_MAGUS/magus.py -d outputs -i $file -o $file.magus
      rm $file
      mv $file.magus $file
      rm -rf outputs
    done
  mv * $DIR_OUTPUT
  cd ..
  rm -rf temp_magus
else 
  cp -r $DIR_INPUT temp1
  mkdir temp2 
  ls temp1/ > loci.list

  cd temp1
  ln -s $DIR_MAFFT/mafft .

  if [ "$MAFFT_METHOD" == "1" ]; then
    ALIGN_fun() {
      ./mafft --thread 1 $1 > ../temp2/$1
    }
  elif [ "$MAFFT_METHOD" == "2" ]; then
    ALIGN_fun() {
      ./mafft --maxiterate 1000 --localpair --thread 1 $1 > ../temp2/$1
    }
  elif [ "$MAFFT_METHOD" == "3" ]; then
    ALIGN_fun() {
      ./mafft --maxiterate 1000 --genafpair --thread 1 $1 > ../temp2/$1
    }
  elif [ "$MAFFT_METHOD" == "4" ]; then
    ALIGN_fun() {
      ./mafft --maxiterate 1000 --globalpair --thread 1 $1 > ../temp2/$1
    }
  else
    echo
  fi

  export -f ALIGN_fun
  cat ../loci.list | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 1 ALIGN_fun %
  cd ..
  mv temp2/* $DIR_OUTPUT/
  rm -rf loci.list temp1 temp2

  for file in $TEMP
    do
      test -s $DIR_OUTPUT/$file && echo "loci $file has been aligned" || $DIR_MAFFT/mafft --thread $THREADS $DIR_INPUT/$file > $DIR_OUTPUT/$file
    done
fi

#Align the the corresponding nucleotide sequences
if [ "$ALIGN_TYPE" == "1" ]; then
echo
read -p "Do you want to align the the corresponding nucleotide sequences? 1. Yes; 2. No:      " NUC_ALIGN
  until [ $NUC_ALIGN -gt 0 -a $NUC_ALIGN -lt 3 ]
    do
      read -p "Do you want to align the the corresponding nucleotide sequences? 1. Yes; 2. No:      " NUC_ALIGN
    done

if [ "$NUC_ALIGN" == "1" ]; then
  mkdir faa
  mv $DIR_OUTPUT/* faa
  mv faa $DIR_OUTPUT/
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
  read -p "Please input the folder names (with its path) containing unaligned nucleotide sequences, for example '2-loci_filter/fna':      " INPUT_TEMP
  INPUT_NUC=$(realpath $(echo $INPUT_TEMP | sed "s/'//g"))
  cd $DIR_OUTPUT/
  rename .faa .fas faa/*.faa
  mkdir temp fna
  cp -r $INPUT_NUC/* temp/
  rename .fna .fas temp/*.fna
  ls faa/ > loci.list
  cat loci.list | $DIR_PARALLEL/parallel -j $THREADS $DIR_TRIMAL/trimal -in faa/{} -out fna/{} -backtrans temp/{}
  rm -rf loci.list temp/
fi
  elif [ "$ALIGN_TYPE" == "2" ]; then
  rename .fna .fas $DIR_OUTPUT/*.fna
fi