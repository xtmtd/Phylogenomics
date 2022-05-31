#!/bin/bash
#2021.08.04 by ZF

#filter loci using gene tree-based methods, i.e. average bootstrap support (ABS), Degree of violation of the molecular clock (DVMC), treeness, signal-to-noise ratio (treeness over rcv), and spurious homologs (possible paralogs, incorrectly assembled sequences) 
#Type 'bash loci_filtering_tree-based.sh'
#parallel, PhyKIT and Seqkit may be required


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


DIR_CURR=$(echo $PWD)

#input the name of input directory
read -p "Please input the name of input directory containing all gene trees, e.g. gene_trees:      " DIR_INPUT_TREE_TEMP
DIR_INPUT_TREE=$(realpath $(echo $DIR_INPUT_TREE_TEMP | sed "s/'//g"))


#input the name of input directory
read -p "Please input the name of input directory containing all loci alignments, e.g. 4-trim/clipkit-kpic:      " DIR_INPUT_ALIGN_TEMP
DIR_INPUT_ALIGN=$(realpath $(echo $DIR_INPUT_ALIGN_TEMP | sed "s/'//g"))


#input the name of output directory
read -p "Please input the name of output directory, or an existing directory:      " DIR_OUTPUT_TEMP
DIR_OUTPUT_TEMP1=$(echo $DIR_OUTPUT_TEMP | sed "s/'//g")
test -d $DIR_OUTPUT_TEMP1 || mkdir -p $DIR_OUTPUT_TEMP1
cd $DIR_OUTPUT_TEMP1 && DIR_OUTPUT=$(echo $PWD) && cd $DIR_CURR


#Check the loci filtering method can be used
read -p "Please input the option for tree-based strategy for loci filtering: 1. average bootstrap support (ABS); 2. degree of violation of the molecular clock (DVMC); 3. treeness; 4. signal-to-noise ratio (treeness over RCV); 5. spurious homologs identification     " FILTER_METHOD
  until [ $FILTER_METHOD -gt 0 -a $FILTER_METHOD -lt 6 ]
    do
      read -p "Please input the option for tree-based strategy for loci filtering: 1. average bootstrap support (ABS); 2. degree of violation of the molecular clock (DVMC); 3. treeness; 4. signal-to-noise ratio (treeness over RCV); 5. spurious homologs identification     " FILTER_METHOD
    done


if [ "$FILTER_METHOD" == "1" ]; then
  #calculate ABS values for trees
  mkdir -p $DIR_OUTPUT/ABS_trees $DIR_OUTPUT/ABS_loci

  cp $DIR_INPUT_TREE/* $DIR_OUTPUT/ABS_trees/
  ls $DIR_INPUT_TREE > $DIR_OUTPUT/loci.list
  cd $DIR_OUTPUT/ABS_trees
  echo -e '\n'
  echo "Calculating average boostrap support (ABS) ......"

  ABS_fun() {
    abs=$(cat $1 | sed "s/)/\n/g" | cut -d ":" -f 1 | sed '1d' | sed "s/;//g" | awk NF | awk '{a+=$1}END{print a/NR}')
    echo -e $1"\t"$abs > $1.abs
  }
  export -f ABS_fun
  cat ../loci.list | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 1 ABS_fun %
  cat *.abs | sort -k2nr > ../ABS.txt
  rm *.abs

  #input the ABS threshod 
  echo -e '\n'
  echo "Read the file ABS.txt in the output folder to determine the minimum ABS threshold (usually 50~90). Loci of higher ABS values are thought to be of more phylogenetic signal."
  echo -e '\n'
  read -p "Please input the ABS threshod:      " THRESHOLD_ABS

  cd $DIR_OUTPUT
  TOTAL_LINE=$(cat ABS.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p ABS.txt | cut -f1)
      num=$(sed -n "$line"p ABS.txt | cut -f2)
      diff=$(echo "scale=4;($num-$THRESHOLD_ABS)"|bc)
      num1=`echo "$diff < 0" |bc`
      test "$num1" = 0 && echo $loci | cut -d "." -f1 >> loci.ABS
    done

  cat loci.list | grep -v -f loci.ABS > loci.ABS.delete
  cat loci.ABS.delete | $DIR_PARALLEL/parallel -j $THREADS rm ABS_trees/{}

  cat loci.ABS | cut -d "." -f1 > temp.list
  cat temp.list | $DIR_PARALLEL/parallel -j $THREADS cp $DIR_INPUT_ALIGN/{}.* ABS_loci/
  rm loci.list loci.ABS.delete temp.list

  echo -e '\n'
  echo "The file loci.ABS contains all the file list of remaining gene trees"
  echo "All the kept gene trees and alignments are deposited in the folder "$DIR_OUTPUT/ABS_trees" and "$DIR_OUTPUT/ABS_loci", respectively"

elif [ "$FILTER_METHOD" == "2" ]; then
  #calculate DVMC values for trees

  #input the file name contaning outgroup taxa names
  read -p "Please input the file name (with its complete path) contaning outgroup taxa names, one name per row, e.g. /home/zf/Desktop/materials/test/6-gene_trees/DVMC/root.txt:      " INPUT_ROOT_TEMP
  INPUT_ROOT=$(realpath $(echo $INPUT_ROOT_TEMP | sed "s/'//g"))

  mkdir -p $DIR_OUTPUT/DVMC_trees $DIR_OUTPUT/DVMC_loci

  cp $DIR_INPUT_TREE/* $DIR_OUTPUT/DVMC_trees/
  ls $DIR_INPUT_TREE > $DIR_OUTPUT/loci.list
  cd $DIR_OUTPUT/DVMC_trees
  echo -e '\n'
  echo "Calculating degree of violation of the molecular clock (DVMC) ......"

  cp $INPUT_ROOT dvmc_root.txt
  ln -s $DIR_PHYKIT/phykit .  
  DVMC_fun() {
    dvmc=$(./phykit dvmc -t $1 -r dvmc_root.txt)
    echo -e $1"\t"$dvmc > $1.dvmc
  }
  export -f DVMC_fun
  cat ../loci.list | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 1 DVMC_fun %
  #for loci in $(cat ../loci.list); do dvmc=$($DIR_PHYKIT/phykit dvmc -t $loci -r $INPUT_ROOT); echo -e $loci"\t"$dvmc > $loci.dvmc; done
  cat *.dvmc | sort -k2n > ../DVMC.txt
  rm *.dvmc dvmc_root.txt phykit

  #input the DVMC threshod 
  echo -e '\n'
  echo "Read the file DVMC.txt in the output folder to determine the maximum DVMC threshold. Lower DVMC values are thought to be desirable because they are indicative of a lower degree of violation in the molecular clock assumption."
  echo -e '\n'
  read -p "Please input the DVMC threshod:      " THRESHOLD_DVMC

  cd $DIR_OUTPUT
  TOTAL_LINE=$(cat DVMC.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p DVMC.txt | cut -f1)
      num=$(sed -n "$line"p DVMC.txt | cut -f2)
      diff=$(echo "scale=4;($num-$THRESHOLD_DVMC)"|bc)
      num1=`echo "$diff >= 0" |bc`
      test "$num1" = 0 && echo $loci >> loci.DVMC
    done

  cat loci.list | grep -v -f loci.DVMC > loci.DVMC.delete
  cat loci.DVMC.delete | $DIR_PARALLEL/parallel -j $THREADS rm DVMC_trees/{}

  cat loci.DVMC | cut -d "." -f1 > temp.list
  cat temp.list | $DIR_PARALLEL/parallel -j $THREADS cp $DIR_INPUT_ALIGN/{}.* DVMC_loci/
  rm loci.list loci.DVMC.delete temp.list

  echo -e '\n'
  echo "The file loci.DVMC contains all the file list of remaining gene trees"
  echo "All the kept gene trees and alignments are deposited in the folder "$DIR_OUTPUT/DVMC_trees" and "$DIR_OUTPUT/DVMC_loci", respectively"

elif [ "$FILTER_METHOD" == "5" ]; then
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

  #remove spurious homolog/outliers
  ls $DIR_INPUT_TREE/ > $DIR_OUTPUT/list.trees
  cp $DIR_INPUT_TREE/* $DIR_OUTPUT/
  cp $DIR_INPUT_ALIGN/* $DIR_OUTPUT/

  cd $DIR_OUTPUT
  mkdir outlier clean_alignments

  #check the outlier sequences/taxa
  ln -s $DIR_PHYKIT/phykit .
  OUTLIER_fun() {
    ./phykit ss $1 -f 20 > $1.outlier
    if [ "$(cat $1.outlier)" = "None" ]; then
      rm $1 $1.outlier
    else
      mv $1.outlier outlier/
      rm $1
    fi
  }
  export -f OUTLIER_fun
  cat list.trees | $DIR_PARALLEL/parallel -j $THREADS OUTLIER_fun
  rm ./phykit

  #remove outlier sequences
  ls outlier/ | cut -d "." -f 1 > list.outlier
  cat list.trees | cut -d "." -f1 | grep -v -f list.outlier > list.kept
  cat list.kept | $DIR_PARALLEL/parallel -j $THREADS mv {}.* clean_alignments/{}.fas
  for loci in $(cat list.outlier)
    do
      cat outlier/$loci*outlier | cut -f 1 > temp
      cat $loci* | $DIR_SEQKIT/seqkit grep -v -f temp > clean_alignments/$loci.fas
      rm $loci*
    done
  rm temp

  echo -e '\n'
  echo "Clean alignments have been deposited in the folder $DIR_OUTPUT/clean_alignments/"
  echo "List of loci containing outlier sequences was saved in the file $DIR_OUTPUT/list.outlier"

else
  #calculate treeness and treeness/RCV values
  mkdir -p $DIR_OUTPUT/tor_trees $DIR_OUTPUT/tor_loci

  cp $DIR_INPUT_TREE/* $DIR_OUTPUT/tor_trees/
  cp $DIR_INPUT_ALIGN/* $DIR_OUTPUT/tor_loci/
  ls $DIR_INPUT_TREE > $DIR_OUTPUT/tree.list
  ls $DIR_INPUT_ALIGN > $DIR_OUTPUT/loci.list

  cp $DIR_INPUT_ALIGN/* $DIR_OUTPUT/tor_trees/
  ls $DIR_OUTPUT/tor_trees > $DIR_OUTPUT/temp.list
  cd $DIR_OUTPUT/tor_trees

  echo -e '\n'
  echo "Calculating treeness/RCV ......"

  ln -s $DIR_PHYKIT/phykit .
  TOR_fun() {
    tor=$(./phykit tor -a $1 -t $2)
    echo -e $1"\t"$tor > $1.tor
  }
  export -f TOR_fun
  cat ../temp.list | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 2 TOR_fun %
  cat *.tor | sort -k3nr > ../TOR.txt
  rm *.tor phykit
  cat ../loci.list | $DIR_PARALLEL/parallel -j $THREADS rm {}

  #input the treeness or treeness/RCV threshod 
  echo -e '\n'
  echo "Read the file TOR.txt in the output folder to determine the minimum threshold. The values of 2nd~4th columns represent treeness/RCV, treeness, and RCV, prespectively. Higher treeness/RCV or treeness values are thought to be desirable."
  echo -e '\n'
  read -p "Please input the threshod:      " THRESHOLD_TOR

  cd $DIR_OUTPUT
  TOTAL_LINE=$(cat TOR.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p TOR.txt | cut -f1)
      if [ "$FILTER_METHOD" == "3" ]; then
        num=$(sed -n "$line"p TOR.txt | awk '{print $3}')
      elif  [ "$FILTER_METHOD" == "4" ]; then
        num=$(sed -n "$line"p TOR.txt | awk '{print $2}')
      fi
      diff=$(echo "scale=4;($num-$THRESHOLD_TOR)"|bc)
      num1=`echo "$diff < 0" |bc`
      test "$num1" = 0 && echo $loci >> loci.TOR
    done

  cat loci.list | grep -v -f loci.TOR > loci.TOR.delete
  cat loci.TOR.delete | $DIR_PARALLEL/parallel -j $THREADS rm tor_loci/{} tor_trees/{}*

  rm loci.list loci.TOR.delete temp.list

  echo -e '\n'
  echo "The file loci.TOR contains all the file list of remaining loci alignments"
  echo "All the kept gene trees and alignments are deposited in the folder "$DIR_OUTPUT/tor_trees and "$DIR_OUTPUT/tor_loci, respectively"

fi

