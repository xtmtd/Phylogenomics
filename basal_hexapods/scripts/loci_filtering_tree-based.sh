#!/bin/bash
#2021.08.04 wrote by ZF
#2023.07.13 revised by DSY
#2023.08.25 by DSY

#filter loci using gene tree-based methods, i.e. average bootstrap support (ABS), Degree of violation of the molecular clock (DVMC), treeness, signal-to-noise ratio (treeness over rcv), spurious homologs (possible paralogs, incorrectly assembled sequences), treeshrink, and inconsistent genes
#Type 'bash loci_filtering_tree-based.sh'
#parallel, R, python, cscvtk, TreShrink, PhyKIT and Seqkit may be required


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


#check csvtk
if [ $(which csvtk) ]
    then
      echo "csvtk ...... OK"
      EXE_CSVTK=$(which csvtk)
      DIR_CSVTK=${EXE_CSVTK%/*}
    else
      until [ -x $DIR_CSVTK/csvtk ]
  do
    read -p "CSVTK is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_CSVTK_TEMP
    DIR_CSVTK=$(realpath $(echo $DIR_CSVTK_TEMP | sed "s/'//g"))
  done
      echo "csvtk ...... OK"
fi


#Check the threads can be used
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
  until [ $THREADS -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done

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
read -p "Please input the option for tree-based strategy for loci filtering: 1. average bootstrap support (ABS); 2. degree of violation of the molecular clock (DVMC); 3. treeness; 4. signal-to-noise ratio (treeness over RCV); 5. spurious homologs identification; 6. treeshrink; 7. inconsistent genes; 8. evolutionary rate; 9. saturation     " FILTER_METHOD
  until [ $FILTER_METHOD -gt 0 -a $FILTER_METHOD -lt 10 ]
    do
      read -p "Please input the option for tree-based strategy for loci filtering: 1. average bootstrap support (ABS); 2. degree of violation of the molecular clock (DVMC); 3. treeness; 4. signal-to-noise ratio (treeness over RCV); 5. spurious homologs identification; 6. treeshrink; 7. inconsistent genes; 8. evolutionary rate; 9. saturation     " FILTER_METHOD
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

  awk '{print $2}' ../ABS.txt > ../temp2
  cat -b ../temp2 | sed "s/ //g" > ../temp
  sed -i '1i\Loci_number\tABS' ../temp
  csvtk -t plot line ../temp -x Loci_number -y ABS --format pdf > ../ABS.pdf
  rm ../temp ../temp2

  #input the ABS threshold 
  echo -e '\n'
  echo "Read the file ABS.txt in the output folder to determine the minimum ABS threshold (usually 50~90). Loci of higher ABS values are thought to be of more phylogenetic signal."
  echo -e '\n'
  read -p "Please input the ABS threshold:      " THRESHOLD_ABS

  cd $DIR_OUTPUT
  TOTAL_LINE=$(cat ABS.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p ABS.txt | cut -f1)
      num=$(sed -n "$line"p ABS.txt | cut -f2)
      diff=$(echo "scale=4;($num-$THRESHOLD_ABS)"|bc)
      num1=`echo "$diff < 0" |bc`
      test "$num1" = 0 && echo $loci >> loci.ABS
    done
 
  cat loci.ABS | sed 's/[ ]//g' > loci.ABS1
  cat loci.ABS1 | sed 's/[\t]//g' > loci.ABS2
  cat loci.list | sed 's/[ ]//g' > loci.list1
  cat loci.list1 | sed 's/[\t]//g' > loci.list2
  rm -rf loci.ABS loci.list loci.ABS1 loci.list1
  mv loci.ABS2 loci.ABS
  mv loci.list2 loci.list
  cp loci.list loci.ABS.delete
  for id in $(cat loci.ABS)
  do
    sed -i "/^$id/d" loci.ABS.delete
  done
  cat loci.ABS.delete | $DIR_PARALLEL/parallel -j $THREADS rm ABS_trees/{}

  sed 's/[\.treefile]\+$//' loci.ABS > temp.list
  cat temp.list | $DIR_PARALLEL/parallel -j $THREADS cp $DIR_INPUT_ALIGN/{} ABS_loci/
  rm loci.list loci.ABS.delete loci.ABS
  mv temp.list loci.ABS

  cd $DIR_OUTPUT/ABS_loci
  $DIR_PHYKIT/phykit cc -a ../loci.ABS -p ABS
  rm *.occupancy
  mv ABS* ..
  AVE=$(echo "scale=2;$(cat $DIR_OUTPUT/ABS.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2)/$(cat $DIR_OUTPUT/loci.ABS | wc -l)"|bc)

  echo -e '\n'
  echo "ABS with $THRESHOLD_ABS (threshold) has preserved $(cat $DIR_OUTPUT/ABS.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/loci.ABS | wc -l) loci (the average site of each locus is $AVE)." | tee -a $DIR_OUTPUT/summary.statistics
  rm ../ABS.partition ../ABS.fa
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

  awk '{print $2}' ../DVMC.txt > ../temp2
  cat -b ../temp2 | sed "s/ //g" > ../temp
  sed -i '1i\Loci_number\tDVMC' ../temp
  csvtk -t plot line ../temp -x Loci_number -y DVMC --format pdf > ../DVMC.pdf
  rm ../temp ../temp2

  #input the DVMC threshold 
  echo -e '\n'
  echo "Read the file DVMC.txt in the output folder to determine the maximum DVMC threshold. Lower DVMC values are thought to be desirable because they are indicative of a lower degree of violation in the molecular clock assumption."
  echo -e '\n'
  read -p "Please input the DVMC threshold:      " THRESHOLD_DVMC

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

  cat loci.DVMC | sed 's/[ ]//g' > loci.DVMC1
  cat loci.DVMC1 | sed 's/[\t]//g' > loci.DVMC2
  cat loci.list | sed 's/[ ]//g' > loci.list1
  cat loci.list1 | sed 's/[\t]//g' > loci.list2
  rm -rf loci.DVMC loci.list loci.DVMC1 loci.list1
  mv loci.DVMC2 loci.DVMC
  mv loci.list2 loci.list
  cp loci.list loci.DVMC.delete
  for id in $(cat loci.DVMC)
  do
    sed -i "/^$id/d" loci.DVMC.delete
  done
  cat loci.DVMC.delete | $DIR_PARALLEL/parallel -j $THREADS rm DVMC_trees/{}

  sed 's/[\.treefile]\+$//' loci.DVMC > temp.list
  cat temp.list | $DIR_PARALLEL/parallel -j $THREADS cp $DIR_INPUT_ALIGN/{} DVMC_loci/
  rm loci.list loci.DVMC.delete loci.DVMC
  mv temp.list loci.DVMC

  cd $DIR_OUTPUT/DVMC_loci
  $DIR_PHYKIT/phykit cc -a ../loci.DVMC -p DVMC
  rm *.occupancy
  mv DVMC* ..
  AVE=$(echo "scale=2;$(cat $DIR_OUTPUT/DVMC.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2)/$(cat $DIR_OUTPUT/loci.DVMC | wc -l)"|bc)

  echo -e '\n'
  echo "DVMC with $THRESHOLD_DVMC (threshold) has preserved $(cat $DIR_OUTPUT/DVMC.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/loci.DVMC | wc -l) loci (the average site of each locus is $AVE)." | tee -a $DIR_OUTPUT/summary.statistics
  rm ../DVMC.partition ../DVMC.fa
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
  rm temp list.trees

  cd $DIR_OUTPUT
  cat list.kept list.outlier > temp
  sed "s/$/&.fas/g" temp > list
  
  for loci in $(cat list)
   do 
     a=$(cat clean_alignments/$loci | wc -l) 
     if [ "$a" == "0" ] 
     then
       echo $loci >> missing.loci
     fi
   done
   sed "s/$/&.fas/g" list.outlier > list.outlier1
   for missing in $(cat missing.loci)
   do
     rm clean_alignments/$missing
     rm outlier/$missing.*
     cat list.outlier1 | sed "s/$missing//g" >> list.outlier2
   done
  
  sort -u list.outlier2 > list.outlier
  ls clean_alignments/ > list
  ls clean_alignments/ | sed "s/.fas/.fas.treefile/g" > list.trees
  cd $DIR_OUTPUT/clean_alignments
  $DIR_PHYKIT/phykit cc -a ../list -p spurious_homolog
  rm *.occupancy
  mv spurious_homolog* ..
  AVE=$(echo "scale=2;$(cat $DIR_OUTPUT/spurious_homolog.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2)/$(cat $DIR_OUTPUT/list | wc -l)"|bc)

  echo -e '\n'
  echo "Spurious_homolog has preserved $(cat $DIR_OUTPUT/spurious_homolog.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/list | wc -l) loci (the average site of each locus is $AVE)." | tee -a $DIR_OUTPUT/summary.statistics
  sed "s/$/&.fas/g" ../list.kept > ../list.kept1
  rm ../spurious_homolog.partition ../spurious_homolog.fa ../list.outlier* ../temp ../list.kept
  mv ../list.kept1 ../list.kept
  cat ../list | grep -v -f ../list.kept > ../list.outlier
  rm ../list
  echo "Clean alignments have been deposited in the folder $DIR_OUTPUT/clean_alignments/"
  echo "List of loci containing outlier sequences was saved in the file $DIR_OUTPUT/list.outlier"

elif [ "$FILTER_METHOD" == "8" ]; then
  #calculate the evolutionary rate
  mkdir -p $DIR_OUTPUT/Rate_trees $DIR_OUTPUT/Rate_loci

  cp $DIR_INPUT_TREE/* $DIR_OUTPUT/Rate_trees/
  ls $DIR_INPUT_TREE > $DIR_OUTPUT/loci.list
  cd $DIR_OUTPUT/Rate_trees
  echo -e '\n'
  echo "Calculating the evolutionary rate ......"

  ln -s $DIR_PHYKIT/phykit .  
  RATE_fun() {
    rate=$(./phykit evo_rate $1)
    echo -e $1"\t"$rate > $1.rate
  }
  export -f RATE_fun
  cat ../loci.list | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 1 RATE_fun %
  cat *.rate | sort -k2n > ../Rate.txt
  rm *.rate phykit

  awk '{print $2}' ../Rate.txt > ../temp2
  cat -b ../temp2 | sed "s/ //g" > ../temp
  sed -i '1i\Loci_number\tEvolutionary_rate' ../temp
  csvtk -t plot line ../temp -x Loci_number -y Evolutionary_rate --format pdf > ../Evolutionary_rate.pdf
  rm ../temp ../temp2

  #input the evolutionary rate threshold 
  echo -e '\n'
  echo "Read the file evolutionary_rate.txt in the output folder to determine the rate threshold. Very high values indicate that those 'fast' genes may be more suitable for shallow phylogenies."
  echo -e '\n'
  read -p "Please input the evolutionary rate threshold:      " THRESHOLD_RATE

  cd $DIR_OUTPUT
  TOTAL_LINE=$(cat Rate.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p Rate.txt | cut -f1)
      num=$(sed -n "$line"p Rate.txt | cut -f2)
      diff=$(echo "scale=4;($num-$THRESHOLD_RATE)"|bc)
      num1=`echo "$diff >= 0" |bc`
      test "$num1" = 0 && echo $loci >> loci.Rate
    done

  cat loci.Rate | sed 's/[ ]//g' > loci.Rate1
  cat loci.Rate1 | sed 's/[\t]//g' > loci.Rate2
  cat loci.list | sed 's/[ ]//g' > loci.list1
  cat loci.list1 | sed 's/[\t]//g' > loci.list2
  rm -rf loci.Rate loci.list loci.Rate1 loci.list1
  mv loci.Rate2 loci.Rate
  mv loci.list2 loci.list
  cp loci.list loci.Rate.delete
  for id in $(cat loci.Rate)
  do
    sed -i "/^$id/d" loci.Rate.delete
  done
  cat loci.Rate.delete | $DIR_PARALLEL/parallel -j $THREADS rm Rate_trees/{}

  sed 's/[\.treefile]\+$//' loci.Rate > temp.list
  cat temp.list | $DIR_PARALLEL/parallel -j $THREADS cp $DIR_INPUT_ALIGN/{} Rate_loci/
  rm loci.list loci.Rate.delete loci.Rate
  mv temp.list loci.Rate

  cd $DIR_OUTPUT/Rate_loci
  $DIR_PHYKIT/phykit cc -a ../loci.Rate -p Rate
  rm *.occupancy
  mv Rate* ..
  AVE=$(echo "scale=2;$(cat $DIR_OUTPUT/Rate.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2)/$(cat $DIR_OUTPUT/loci.Rate | wc -l)"|bc)

  echo -e '\n'
  echo "Evolutionary rate with $THRESHOLD_RATE (threshold) has preserved $(cat $DIR_OUTPUT/Rate.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/loci.Rate | wc -l) loci (the average site of each locus is $AVE)." | tee -a $DIR_OUTPUT/summary.statistics
  rm ../Rate.partition ../Rate.fa
  echo "The file loci.Rate contains all the file list of remaining gene trees"
  echo "All the kept gene trees and alignments are deposited in the folder "$DIR_OUTPUT/Rate_trees" and "$DIR_OUTPUT/Rate_loci", respectively"

elif [ "$FILTER_METHOD" == "6" ]; then
#check treeshrink
if [ $(which run_treeshrink.py) ]
    then
      echo "run_treeshrink.py ...... OK"
      EXE_TREESHRINK=$(which run_treeshrink.py)
      DIR_TREESHRINK=${EXE_TREESHRINK%/*}
    else
      until [ -x $DIR_TREESHRINK/run_treeshrink.py ]
        do
          read -p "run_treeshrink.py is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/TreeShrink-1.3.8b):      " DIR_TREESHRINK_TEMP
          DIR_TREESHRINK=$(realpath $(echo $DIR_TREESHRINK_TEMP | sed "s/'//g"))
        done
      echo "treeshrink ...... OK"
fi

  #input the number of α threshold (α = 0.05 is default)
  read -p "Please input the number of α threshold. Default is 0.05 (e.g. 0.05):      " THRESHOLD
    until [ $THRESHOLD ]
      do
        read -p "Please input the number of α threshold. Default is 0.05 (e.g. 0.05):      " THRESHOLD
      done

  #copy the alignments and gene trees, and extract the loci
  ls $DIR_INPUT_ALIGN/ > $DIR_OUTPUT/loci.list
  mkdir -p $DIR_OUTPUT/outlier && cd $DIR_OUTPUT/outlier
  cat ../loci.list | $DIR_PARALLEL/parallel -j $THREADS mkdir {}
  cat ../loci.list | $DIR_PARALLEL/parallel -j $THREADS cp $DIR_INPUT_ALIGN/{} {}
  cat ../loci.list | $DIR_PARALLEL/parallel -j $THREADS cp $DIR_INPUT_TREE/{}.treefile {}
  cd ..

  #treeshrink analysis
  TREESHRINK_fun() {
    $DIR_TREESHRINK/run_treeshrink.py -i outlier/ -t $1.treefile -q "$THRESHOLD" -a $1 -O $1.treeshrink
    }
  export -f TREESHRINK_fun
  cat loci.list | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 1 TREESHRINK_fun %
  
  #extract filtered sequences
  for loci in $(cat loci.list)
  do
    if test -s outlier/$loci/*.txt; then
      echo "$loci" >> list.outlier
    else 
      echo "$loci" >> list.keep
    fi
  done

  mkdir temp temp_trees
  cp $DIR_INPUT_TREE/* temp_trees
  PRUNE_fun() {
    $DIR_PHYKIT/phykit prune_tree temp_trees/$1.treefile outlier/$1/$1.treeshrink.txt -o temp/$1.treeshrink.treefile
   }
  export -f PRUNE_fun
  cat list.outlier | $DIR_PARALLEL/parallel -j $THREADS PRUNE_fun
 
  mkdir treeshrink_loci treeshrink_trees
  cat list.keep | $DIR_PARALLEL/parallel -j $THREADS cp outlier/{}/{} treeshrink_loci/{}
  cat list.keep | $DIR_PARALLEL/parallel -j $THREADS cp $DIR_INPUT_TREE/{}.treefile treeshrink_trees/{}
  
  cat list.outlier | wc -l > temp1
  if test -s temp1; then
    cat list.outlier | $DIR_PARALLEL/parallel -j $THREADS cp outlier/{}/{} treeshrink_loci/{}
    cat list.outlier | $DIR_PARALLEL/parallel -j $THREADS mv temp/{}.treeshrink.treefile treeshrink_trees/{}

  rename 's/.fas/.fas.treefile/g' treeshrink_trees/*

  cd $DIR_OUTPUT/treeshrink_loci
  $DIR_PHYKIT/phykit cc -a ../loci.list -p treeshrink
  rm *.occupancy
  mv treeshrink* ..
  AVE=$(echo "scale=2;$(cat $DIR_OUTPUT/treeshrink.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2)/$(cat $DIR_OUTPUT/loci.list | wc -l)"|bc)
  cd ..

  echo -e '\n'
  echo "Treeshrink $THRESHOLD (α threshold) has preserved $(cat $DIR_OUTPUT/treeshrink.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/loci.list | wc -l) loci (the average site of each locus is $AVE)." | tee -a $DIR_OUTPUT/summary.statistics
  rm treeshrink.partition treeshrink.fa
  echo "List of loci containing outlier sequences was saved in the file $DIR_OUTPUT/list.outlier."
  echo "All alignments and gene trees have been deposited in the folder $DIR_OUTPUT/treeshrink_loci and $DIR_OUTPUT/treeshrink_trees, respectively."
  echo "The taxa of each alignment must more than 10. If not, these alignments will not be analysed. Check them carefully!"
  else
  rm -rf list.outlier
  echo -e '\n'
  echo "Congratulations! All alignments have no abnormally long branches."
  echo -e '\n'
  echo "Treeshrink $THRESHOLD (α threshold) has preserved $(cat $DIR_OUTPUT/treeshrink.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/loci.list | wc -l) loci (the average site of each locus is $AVE)." | tee -a $DIR_OUTPUT/summary.statistics
  rm treeshrink.partition treeshrink.fa
  echo "All alignments and gene trees have been deposited in the folder $DIR_OUTPUT/treeshrink_loci and $DIR_OUTPUT/treeshrink_trees, respectively."
  echo "The taxa of each alignment must more than 10. If not, these alignments will not be analysed. Check them carefully!"
  fi
  rm -rf outlier/ temp1 temp/ temp_trees/


elif [ "$FILTER_METHOD" == "7" ]; then
  #check ASTRAL
  until [ -s $DIR_ASTRAL/astral*.jar ]
     do
       read -p "ASTRAL executable file is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/Astral-5.7.8):      " DIR_ASTRAL_TEMP
       DIR_ASTRAL=$(realpath $(echo $DIR_ASTRAL_TEMP | sed "s/'//g"))
     done
  echo "ASTRAL ...... OK"

  #input the resulting topological in supermatrix-based (T1) phylogeny
  echo
  read -p "Please input the resulting topological in supermatrix-based (T1) phylogeny:    " T1_TEMP
  T1=$(realpath $(echo $T1_TEMP | sed "s/'//g"))
  
  #input the resulting topological in supertree-based (T2) phylogeny
  echo
  read -p "Please input the resulting topological in supertree-based (T2) phylogeny:    " T2_TEMP
  T2=$(realpath $(echo $T2_TEMP | sed "s/'//g"))

  mkdir -p $DIR_OUTPUT/GLS $DIR_OUTPUT/GQS 
  #Estimate site-wise loglk for alternative hypothesis
  cd $DIR_OUTPUT/GLS
  mkdir fas && cd fas
  cp $DIR_INPUT_ALIGN/* ./
  ls > ../loci.alignments
  $DIR_PHYKIT/phykit cc -a ../loci.alignments -p GLS
  mv GLS.fa ../supermatrix.fa
  mv GLS.partition ..
  cd ..
  rm -r fas/
  cat $T1 $T2 > ML2ASTRAL.tre
  sed 's/.*,/LG,/g' GLS.partition > partition.txt   #set model in parition file as 'LG' instead of "AUTO"
  iqtree -s supermatrix.fa -Q partition.txt -m LG+F+R4 -z ML2ASTRAL.tre -wsl -T $THREADS --prefix ML2ASTRAL
  cat partition.txt | awk '{print $2}' | sed 's/=/\t/g' | sed 's/-/\t/g' > temp1
  cat temp1 | awk '{print $1}' > gene_id
  cat ML2ASTRAL.sitelh | sed -n 2p | awk '{$1=null;print $0}' | sed 's/ /\t/g' | csvtk -t transpose > Tree1
  cat ML2ASTRAL.sitelh | sed -n 3p | awk '{$1=null;print $0}' | sed 's/ /\t/g' | csvtk -t transpose > Tree2
  sed -i '1d' Tree1
  sed -i '1d' Tree2
  for gene in $(cat gene_id)
  do
    i=$(cat temp1 | grep "^$gene" | awk '{print $2}')
    j=$(cat temp1 | grep "^$gene" | awk '{print $3}')
    sed -n "$i, $j p" Tree1 > $gene.tr1
    sed -n "$i, $j p" Tree2 > $gene.tr2
    cat $gene.tr1 | awk '{sum+=$1}END{print sum}' >> tr1_GLS
    cat $gene.tr2 | awk '{sum+=$1}END{print sum}' >> tr2_GLS
  done
  paste tr1_GLS tr2_GLS > GLS
  awk '{$3 = $1 - $2}1' GLS > tr_GLS
  n=$(cat tr_GLS | wc -l)
  for i in `seq 1 1 $n`
  do
    a=$(sed -n "$i p" tr1_GLS)
    b=$(sed -n "$i p" tr2_GLS)
    if [ `expr $a \> $b` -eq 0 ]; then
      echo tr1 >> tree_supported
    else
      echo tr2 >> tree_supported
    fi
  done
  paste gene_id tree_supported tr_GLS > GLS_table
  sed -i '1i gene_id\ttree_supported\ttr1_GLS\ttr2_GLS\tdf_GLS' GLS_table
  cat GLS_table | sed 's/ /\t/g' > GLS_table.txt
  #table format 'gene_id  tree_supported  tr1_GLS tr2_GLS df_GLS'

#estimate the quartet score for alternative hypothesis
cd $DIR_OUTPUT/GQS
mkdir gene_trees
#prune tips of T1 and T2 to make their tips equal to gene trees
ls $DIR_INPUT_TREE/ > tree.list
cp $DIR_INPUT_TREE/* gene_trees/
cp $T1 T1.tre
cp $T2 T2.tre
ln -s $DIR_ASTRAL/*.jar .
ln -s $DIR_PHYKIT/phykit .
GQS_fun() {
  ./phykit tl gene_trees/$1 > tip.$1
  ./phykit tl T1.tre > tip.T1.$1
  cat tip.T1.$1 | grep -v -f tip.$1 > tip.missing.$1
  ./phykit prune T1.tre tip.missing.$1 -o T1.prune.tre.$1
  ./phykit prune T2.tre tip.missing.$1 -o T2.prune.tre.$1
  java -jar ./*.jar -i gene_trees/$1 -q T1.prune.tre.$1 2> log.txt.T1.$1
  GQS_T1=$(cat log.txt.T1.$1 | grep "^Final quartet" | cut -d ":" -f2 | sed "s/ //g")
  java -jar ./*.jar -i gene_trees/$1 -q T2.prune.tre.$1 2> log.txt.T2.$1
  GQS_T2=$(cat log.txt.T2.$1 | grep "^Final quartet" | cut -d ":" -f2 | sed "s/ //g")
  diff=$(echo "scale=0;($GQS_T1-$GQS_T2)"|bc)
  echo -e $1"\t"$GQS_T1"\t"$GQS_T2"\t"$diff > gene_trees/$1.GQS
  rm tip*$1 log.txt*$1 *prune.tre.$1
}
export -f GQS_fun
cat tree.list | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 1 GQS_fun %
cat gene_trees/*GQS > GQS_table.txt
rm phykit
#table format 'tree_id  tr1_GLS tr2_GLS df_GLS'

#get te inconsistent loci
cd $DIR_OUTPUT/
cat GQS/GQS_table.txt | cut -f1 | cut -d "." -f1 > loci.list
for loci in $(cat loci.list);
  do
    num1=$(cat GLS/GLS_table.txt | grep $loci | cut -f5 | sed "s/e-//g;s/e+//g")
    num2=$(cat GQS/GQS_table.txt | grep $loci | cut -f4)
    c_num1=`echo "$num1 > 0" |bc`
    c_num2=`echo "$num2 > 0" |bc`
    test "$c_num1" = 1 -a "$c_num2" = 1 && echo $loci >> loci.consistent
    c_num3=`echo "$num1 < 0" |bc`
    c_num4=`echo "$num2 < 0" |bc`
    test "$c_num3" = 1 -a "$c_num4" = 1 && echo $loci >> loci.consistent
  done
  mkdir consistent_loci consistent_trees
  for loci in $(cat loci.consistent)
  do
    cp $DIR_INPUT_ALIGN/$loci.* consistent_loci
    cp $DIR_INPUT_TREE/$loci.* consistent_trees
  done
  
  rm loci.consistent
  cd $DIR_OUTPUT/consistent_loci
  ls > ../loci.consistent
  $DIR_PHYKIT/phykit cc -a ../loci.consistent -p consistent
  rm *.occupancy
  mv consistent* ..
  AVE=$(echo "scale=2;$(cat $DIR_OUTPUT/consistent.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2)/$(cat $DIR_OUTPUT/loci.consistent | wc -l)"|bc)

  echo -e '\n'
  echo "Consistent gene has preserved $(cat $DIR_OUTPUT/consistent.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/loci.consistent | wc -l) loci (the average site of each locus is $AVE)." | tee -a $DIR_OUTPUT/summary.statistics
  mv ../GLS/GLS_table.txt ..
  mv ../GQS/GQS_table.txt ..
  rm -rf ../consistent.partition ../consistent.fa ../GLS/ ../GQS/ ../loci.list
  echo "The file loci.consistent contains all the file list of remaining gene trees"
  echo "All the kept gene trees and alignments are deposited in the folder "$DIR_OUTPUT/consistent_trees" and "$DIR_OUTPUT/consistent_loci", respectively"


elif [ "$FILTER_METHOD" == "9" ]; then
  #calculate the saturation
  mkdir -p $DIR_OUTPUT/Saturation_trees $DIR_OUTPUT/Saturation_loci

  cp $DIR_INPUT_TREE/* $DIR_OUTPUT/Saturation_trees/
  cp $DIR_INPUT_ALIGN/* $DIR_OUTPUT/Saturation_loci/
  ls $DIR_INPUT_ALIGN > $DIR_OUTPUT/loci.list

  cp $DIR_INPUT_ALIGN/* $DIR_OUTPUT/Saturation_trees/
  ls $DIR_OUTPUT/Saturation_trees > $DIR_OUTPUT/temp.list
  cd $DIR_OUTPUT/Saturation_trees

  echo -e '\n'
  echo "Calculating the saturation ......"

  ln -s $DIR_PHYKIT/phykit .
  SAT_fun() {
    saturation=$(./phykit sat -a $1 -t $2)
    echo -e $1"\t"$saturation > $1.saturation
  }
  export -f SAT_fun
  cat ../temp.list | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 2 SAT_fun %
  cat *.saturation | sort -k2n > ../saturation.txt
  rm *.saturation phykit
  cat ../loci.list | $DIR_PARALLEL/parallel -j $THREADS rm {}

  awk '{print $2}' ../saturation.txt > ../temp2
  cat -b ../temp2 | sed "s/ //g" > ../temp
  sed -i '1i\Loci_number\tSaturation' ../temp
  csvtk -t plot line ../temp -x Loci_number -y Saturation --format pdf > ../Saturation.pdf
  rm ../temp ../temp2

  #input the value of saturation
  echo -e '\n'
  echo "Read the file saturation.txt in the output folder to determine the minimum the value of saturation. Loci of higher values of saturation are thought to be desirable."
  echo -e '\n'
  read -p "Please input the value of saturation:      " THRESHOLD_SAT

  cd $DIR_OUTPUT
  TOTAL_LINE=$(cat saturation.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p saturation.txt | cut -f1)
      num=$(sed -n "$line"p saturation.txt | cut -f2)
      diff=$(echo "scale=4;($num-$THRESHOLD_SAT)"|bc)
      num1=`echo "$diff < 0" |bc`
      test "$num1" = 0 && echo $loci >> loci.saturation
    done
  
  rm loci.list
  ls $DIR_INPUT_TREE > $DIR_OUTPUT/loci.list

  cat loci.saturation | sed 's/[ ]//g' > loci.saturation1
  cat loci.saturation1 | sed 's/[\t]//g' > loci.saturation2
  cat loci.list | sed 's/[ ]//g' > loci.list1
  cat loci.list1 | sed 's/[\t]//g' > loci.list2
  rm -rf loci.saturation loci.list loci.saturation1 loci.list1
  mv loci.saturation2 loci.saturation
  mv loci.list2 loci.list
  cp loci.list loci.saturation.delete
  for id in $(cat loci.saturation)
  do
    sed -i "/^$id/d" loci.saturation.delete
  done
  cat loci.saturation.delete | $DIR_PARALLEL/parallel -j $THREADS rm Saturation_trees/{}

  sed 's/[\.treefile]\+$//' loci.saturation > temp.list
  cat temp.list | $DIR_PARALLEL/parallel -j $THREADS cp $DIR_INPUT_ALIGN/{} Saturation_loci/
  rm loci.list loci.saturation.delete loci.saturation
  mv temp.list loci.saturation

  cd $DIR_OUTPUT/Saturation_loci
  $DIR_PHYKIT/phykit cc -a ../loci.saturation -p saturation
  rm *.occupancy
  mv saturation* ..
  AVE=$(echo "scale=2;$(cat $DIR_OUTPUT/saturation.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2)/$(cat $DIR_OUTPUT/loci.saturation | wc -l)"|bc)

  echo -e '\n'
  echo "Saturation with $THRESHOLD_SAT (threshold) has preserved $(cat $DIR_OUTPUT/saturation.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/loci.saturation | wc -l) loci (the average site of each locus is $AVE)." | tee -a $DIR_OUTPUT/summary.statistics
  rm ../saturation.partition ../saturation.fa
  echo "The file loci.saturation contains all the file list of remaining gene trees"
  echo "All the kept gene trees and alignments are deposited in the folder "$DIR_OUTPUT/Saturation_trees" and "$DIR_OUTPUT/Saturation_loci", respectively"


else
  #calculate treeness and treeness/RCV values
  mkdir -p $DIR_OUTPUT/tor_trees $DIR_OUTPUT/tor_loci

  cp $DIR_INPUT_TREE/* $DIR_OUTPUT/tor_trees/
  cp $DIR_INPUT_ALIGN/* $DIR_OUTPUT/tor_loci/
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

  if [ "$FILTER_METHOD" == "3" ]; then
  awk '{print $3}' ../TOR.txt > ../temp1
  cat ../temp1 | csvtk -t -H sort -k1:n > ../temp2
  cat -b ../temp2 | sed "s/ //g" > ../temp
  sed -i '1i\Loci_number\ttreeness' ../temp
  csvtk -t plot line ../temp -x Loci_number -y treeness --format pdf > ../treeness.pdf
  rm ../temp ../temp2 ../temp1
  else
  awk '{print $2}' ../TOR.txt > ../temp1
  cat ../temp1 | csvtk -t -H sort -k1:n > ../temp2
  cat -b ../temp2 | sed "s/ //g" > ../temp
  sed -i '1i\Loci_number\ttreeness_over_RCV' ../temp
  csvtk -t plot line ../temp -x Loci_number -y treeness_over_RCV --format pdf > ../treeness_over_RCV.pdf
  rm ../temp ../temp2 ../temp1
  fi

  #input the treeness or treeness/RCV threshold 
  echo -e '\n'
  echo "Read the file TOR.txt in the output folder to determine the minimum threshold. The values of 2nd~4th columns represent treeness/RCV, treeness, and RCV, prespectively. Higher treeness/RCV or treeness values are thought to be desirable."
  echo -e '\n'
  read -p "Please input the threshold:      " THRESHOLD_TOR

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
  
  cat loci.TOR | sed 's/[ ]//g' > loci.TOR1
  cat loci.TOR1 | sed 's/[\t]//g' > loci.TOR2
  cat loci.list | sed 's/[ ]//g' > loci.list1
  cat loci.list1 | sed 's/[\t]//g' > loci.list2
  rm -rf loci.TOR loci.list loci.TOR1 loci.list1
  mv loci.TOR2 loci.TOR
  mv loci.list2 loci.list
  cp loci.list loci.TOR.delete
  for id in $(cat loci.TOR)
  do
    sed -i "/^$id/d" loci.TOR.delete
  done
  cat loci.TOR.delete | $DIR_PARALLEL/parallel -j $THREADS rm tor_loci/{} tor_trees/{}*

  rm loci.list loci.TOR.delete temp.list

  cd $DIR_OUTPUT/tor_loci
  $DIR_PHYKIT/phykit cc -a ../loci.TOR -p TOR
  rm *.occupancy
  mv TOR* ..
  AVE=$(echo "scale=2;$(cat $DIR_OUTPUT/TOR.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2)/$(cat $DIR_OUTPUT/loci.TOR | wc -l)"|bc)

  echo -e '\n'
  echo "Treeness with $THRESHOLD_TOR (threshold) has preserved $(cat $DIR_OUTPUT/TOR.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/loci.TOR | wc -l) loci (the average site of each locus is $AVE)." | tee -a $DIR_OUTPUT/summary.statistics
  rm ../TOR.partition ../TOR.fa
  echo "The file loci.TOR contains all the file list of remaining loci alignments"
  echo "All the kept gene trees and alignments are deposited in the folder "$DIR_OUTPUT/tor_trees and "$DIR_OUTPUT/tor_loci, respectively"

fi