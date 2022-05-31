#!/bin/bash
#2021.08.04 by ZF

#filter loci by detecting alignment length, number of parsimony-informative sites, percentage of parsimony-informative sites in the alignment, GC content, compositional heterogeneity (RCV, Relative composition variability), evolutionary rate (average pairwise identity) and likelihood mapping criteria.
#Type 'bash loci_filtering_alignment-based.sh'
#parallel, PhyKIT and IQ-TREE may be required


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
read -p "Please input the name of input directory containing all alignments, e.g. 4-trim/clipkit-kpic:      " DIR_INPUT_TEMP
DIR_INPUT=$(realpath $(echo $DIR_INPUT_TEMP | sed "s/'//g"))


#input the name of output directory
read -p "Please input the name of output directory, or an existing directory:      " DIR_OUTPUT_TEMP
DIR_OUTPUT_TEMP1=$(echo $DIR_OUTPUT_TEMP | sed "s/'//g")
test -d $DIR_OUTPUT_TEMP1 || mkdir -p $DIR_OUTPUT_TEMP1
cd $DIR_OUTPUT_TEMP1 && DIR_OUTPUT=$(echo $PWD) && cd $DIR_CURR


#Check the loci filtering method can be used
read -p "Please input the option for alignment-based strategy for loci filtering: 1. alignment length; 2. number of parsimony-informative sites; 3. percentage of parsimony-informative sites in the alignment; 4. GC content; 5. RCV (Relative composition variability); 6. evolutionary rate (average pairwise identity); 7. likelihood mapping; 8. symmetry tests against SRH hypotheses     " FILTER_METHOD
  until [ $FILTER_METHOD -gt 0 -a $FILTER_METHOD -lt 9 ]
    do
      read -p "Please input the option for alignment-based strategy for loci filtering: 1. alignment length; 2. number of parsimony-informative sites; 3. percentage of parsimony-informative sites in the alignment; 4. GC content; 5. RCV (Relative composition variability); 6. evolutionary rate (average pairwise identity); 7. likelihood mapping; 8. symmetry tests against SRH hypotheses     " FILTER_METHOD
    done



if [ "$FILTER_METHOD" == "5" ]; then
  #calculate RCV values for alignments
  mkdir -p $DIR_OUTPUT/rcv

  cp $DIR_INPUT/* $DIR_OUTPUT/rcv/
  ls $DIR_INPUT > $DIR_OUTPUT/loci.list
  cd $DIR_OUTPUT/rcv
  echo -e '\n'
  echo "Calculating RCV ......"

  ln -s $DIR_PHYKIT/phykit .
  RCV_fun() {
    rcv=$(./phykit rcv $1)
    echo -e $1"\t"$rcv > $1.rcv
  }
  export -f RCV_fun
  cat ../loci.list | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 1 RCV_fun %
  cat *.rcv | sort -k2n > ../rcv.txt
  rm *.rcv phykit

  #input the RCV threshod 
  echo -e '\n'
  echo "Read the file rcv.txt in the output folder to determine the RCV threshold. Lower RCV values are thought to be desirable because they represent a lower composition bias in an alignment."
  echo -e '\n'
  read -p "Please input the RCV threshod:      " THRESHOLD_RCV

  TOTAL_LINE=$(cat ../rcv.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p ../rcv.txt | cut -f1)
      num=$(sed -n "$line"p ../rcv.txt | cut -f2)
      diff=$(echo "scale=4;($num-$THRESHOLD_RCV)"|bc)
      num1=`echo "$diff >= 0" |bc`
      test "$num1" = 0 && echo $loci >> ../loci.rcv
    done

  cat ../loci.list | grep -v -f ../loci.rcv > ../loci.rcv.delete
  cat ../loci.rcv.delete | $DIR_PARALLEL/parallel -j $THREADS rm {}
  cd ..
  rm loci.list loci.rcv.delete

  echo -e '\n'
  echo "The file loci.rcv contains all the file list of remaining alignments."
  echo "All the kept alignments are deposited in the folder 'rcv'."

elif [ "$FILTER_METHOD" == "4" ]; then
  #calculate GC values for nucleotide alignments
  mkdir -p $DIR_OUTPUT/GC

  cp $DIR_INPUT/* $DIR_OUTPUT/GC/
  ls $DIR_INPUT > $DIR_OUTPUT/loci.list
  cd $DIR_OUTPUT/GC
  echo -e '\n'
  echo "Calculating GC content ......"

  ln -s $DIR_PHYKIT/phykit .
  GC_fun() {
    gc=$(./phykit gc $1)
    echo -e $1"\t"$gc > $1.gc
  }
  export -f GC_fun
  cat ../loci.list | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 1 GC_fun %
  cat *.gc | sort -k2n > ../GC.txt
  rm *.gc phykit

  #input the GC threshod 
  echo -e '\n'
  echo "Read the file GC.txt in the output folder to determine the GC percentage threshold. Very high (often > 0.8) GC values are thought to be unreliable."
  echo -e '\n'
  read -p "Please input the GC percentage threshod:      " THRESHOLD_GC

  TOTAL_LINE=$(cat ../GC.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p ../GC.txt | cut -f1)
      num=$(sed -n "$line"p ../GC.txt | cut -f2)
      diff=$(echo "scale=4;($num-$THRESHOLD_GC)"|bc)
      num1=`echo "$diff >= 0" |bc`
      test "$num1" = 0 && echo $loci >> ../loci.GC
    done

  cat ../loci.list | grep -v -f ../loci.GC > ../loci.GC.delete
  cat ../loci.GC.delete | $DIR_PARALLEL/parallel -j $THREADS rm {}
  cd ..
  rm loci.list loci.GC.delete

  echo -e '\n'
  echo "The file loci.GC contains all the file list of remaining alignments"
  echo "All the kept alignments are deposited in the folder 'GC'"

elif [ "$FILTER_METHOD" == "6" ]; then
  #calculate values of mean average pairwise identity
  mkdir -p $DIR_OUTPUT/evolutionary_rate

  cp $DIR_INPUT/* $DIR_OUTPUT/evolutionary_rate/
  ls $DIR_INPUT > $DIR_OUTPUT/loci.list
  cd $DIR_OUTPUT/evolutionary_rate
  echo -e '\n'
  echo "Calculating the mean pairwise identity ......"

  ln -s $DIR_PHYKIT/phykit .
  RATE_fun() {
    rate=$(./phykit pi $1 | head -n 1 | sed "s/mean: //g")
    echo -e $1"\t"$rate > $1.rate
  }
  export -f RATE_fun
  cat ../loci.list | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 1 RATE_fun %
  cat *.rate | sort -k2n > ../evolutionary_rate.txt
  rm *.rate phykit

  #input the evolutionary rate threshod 
  echo -e '\n'
  echo "Read the file evolutionary_rate.txt in the output folder to determine the rate threshold. Very high values indicate that those 'fast' genes may be more suitable for shallow phylogenies."
  echo -e '\n'
  read -p "Please input the evolutionary_rate threshod:      " THRESHOLD_RATE

  TOTAL_LINE=$(cat ../evolutionary_rate.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p ../evolutionary_rate.txt | cut -f1)
      num=$(sed -n "$line"p ../evolutionary_rate.txt | cut -f2)
      diff=$(echo "scale=4;($num-$THRESHOLD_RATE)"|bc)
      num1=`echo "$diff >= 0" |bc`
      test "$num1" = 0 && echo $loci >> ../loci.evolutionary_rate
    done

  cat ../loci.list | grep -v -f ../loci.evolutionary_rate > ../loci.evolutionary_rate.delete
  cat ../loci.evolutionary_rate.delete | $DIR_PARALLEL/parallel -j $THREADS rm {}
  cd ..
  rm loci.list loci.evolutionary_rate.delete

  echo -e '\n'
  echo "The file loci.evolutionary_rate contains all the file list of remaining alignments."
  echo "All the kept alignments are deposited in the folder 'evolutionary_rate.'"


elif [ "$FILTER_METHOD" == "7" ]; then
  #calculate number of fully resolved  quartets for each locus

  #check IQ-TREE
  if [ $(which iqtree) ]
    then
      echo "IQ-TREE ...... OK"
      EXE_IQTREE=$(which iqtree)
      DIR_IQTREE=${EXE_IQTREE%/*}
    else
      until [ -x $EXE_IQTREE ]
        do
          read -p "IQ-TREE is not found. Please input the location of its executive file (e.g. /usr/local/bin/iqtree2):      " EXE_IQTREE_TEMP
          EXE_IQTREE=$(realpath $(echo $EXE_IQTREE_TEMP | sed "s/'//g"))
          DIR_IQTREE=${EXE_IQTREE%/*}
        done
      echo "IQ-TREE ...... OK"
  fi

  #Check the type of input alignments
  read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide   " ALIGN_TYPE
    until [ $ALIGN_TYPE -gt 0 -a $ALIGN_TYPE -lt 3 ]
      do
        read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide      " ALIGN_TYPE
      done

  mkdir -p $DIR_OUTPUT/likelihood_mapping $DIR_OUTPUT/iqtree $DIR_OUTPUT/tempLM

  cp $DIR_INPUT/* $DIR_OUTPUT/tempLM/
  ls $DIR_INPUT > $DIR_OUTPUT/loci.list
  cd $DIR_OUTPUT/tempLM
  echo -e '\n'
  echo "Calculating the number of fully resolved  quartets for each locus using likelihood mapping method ......"

  ln -s $EXE_IQTREE iqtree

#Check the substitution model
if [ "$ALIGN_TYPE" == "1" ]; then
  read -p "Please input the option for protein substitution model: 1. automatically ModelFinder searching; 2. restricted to LG model to reduce computational burden (-mset LG); 3. mixture model EX_EHO (-mset EX_EHO, may be very time-consuming but more accurate); 4. insect universal model Q.insect (may be better than LG model for insects); 5. mitochondrial (--msub mitochondrial)         " MODEL_INPUT
    until [ $MODEL_INPUT -gt 0 -a $MODEL_INPUT -lt 6 ]
      do
        read -p "Please input the option for protein substitution model: 1. automatically ModelFinder searching; 2. restricted to LG model to reduce computational burden (-mset LG); 3. mixture model EX_EHO (-mset EX_EHO, may be very time-consuming but more accurate); 4. insect universal model Q.insect (may be better than LG model for insects); 5. mitochondrial (--msub mitochondrial)         " MODEL_INPUT
      done
  [ "$MODEL_INPUT" == "1" ] && MODEL=""
  [ "$MODEL_INPUT" == "2" ] && MODEL="--mset LG"
  [ "$MODEL_INPUT" == "3" ] && MODEL="--mset EX_EHO"
  [ "$MODEL_INPUT" == "4" ] && MODEL="--mset Q.insect"
  [ "$MODEL_INPUT" == "5" ] && MODEL="--msub mitochondrial"

  cat ../loci.list | $DIR_PARALLEL/parallel -j $THREADS ./iqtree -s {} $MODEL -lmap ALL -n 0

else
  cat ../loci.list | $DIR_PARALLEL/parallel -j $THREADS ./iqtree -s {} -mset GTR -lmap ALL -n 0
fi

  for loci in $(cat ../loci.list)
    do
      quartet=$(cat $loci.iqtree | grep "Number of fully resolved  quartets" | cut -d "=" -f2 | sed "s/%)//g")
      echo -e $loci"\t"$quartet > $loci.quartet
    done
  cat *.quartet | sort -k2nr > ../likelihood_mampping_quartet.txt
  mv *.iqtree ../iqtree/
  cd .. 
  rm -rf tempLM

  #input the quartet threshod 
  echo -e '\n'
  echo "Read the file likelihood_mampping_quartet.txt in the output folder to determine the quartet threshold. High values indicate that strong phylogenetic informativeness/signal."
  echo -e '\n'
  read -p "Please input the quartet threshod:      " THRESHOLD_QUARTET

  TOTAL_LINE=$(cat likelihood_mampping_quartet.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p likelihood_mampping_quartet.txt | cut -f1)
      num=$(sed -n "$line"p likelihood_mampping_quartet.txt | cut -f2)
      diff=$(echo "scale=2;($num-$THRESHOLD_QUARTET)"|bc)
      num1=`echo "$diff >= 0" |bc`
      test "$num1" = 1 && echo $loci >> loci.likelihood_mampping
    done

  cat loci.likelihood_mampping | $DIR_PARALLEL/parallel -j $THREADS cp $DIR_INPUT/{} likelihood_mapping/

  echo -e '\n'
  echo "The file loci.likelihood_mampping contains all the file list of remaining alignments."
  echo "All the kept alignments are deposited in the folder 'likelihood_mampping.'"


elif [ "$FILTER_METHOD" == "8" ]; then
  #symmetry tests

  #check IQ-TREE
  if [ $(which iqtree) ]
    then
      echo "IQ-TREE ...... OK"
      EXE_IQTREE=$(which iqtree)
      DIR_IQTREE=${EXE_IQTREE%/*}
    else
      until [ -x $EXE_IQTREE ]
        do
          read -p "IQ-TREE is not found. Please input the location of its executive file (e.g. /usr/local/bin/iqtree2):      " EXE_IQTREE_TEMP
          EXE_IQTREE=$(realpath $(echo $EXE_IQTREE_TEMP | sed "s/'//g"))
          DIR_IQTREE=${EXE_IQTREE%/*}
        done
      echo "IQ-TREE ...... OK"
  fi

  #generate concatenated matrix
  cd $DIR_OUTPUT
  cp -r $DIR_INPUT symtest
  ls symtest/ > all.alignments
  cd symtest/
  $DIR_PHYKIT/phykit cc -a ../all.alignments -p all
  mv all.* ..
  cd .. 
  rm all.alignments all.occupancy
  
  #calculate symtest
  ln -s $EXE_IQTREE iqtree
  iqtree -s all.fa -p all.partition --symtest-only -T $THREADS
  cat all.partition.symtest.csv | grep -v "^#" | cut -d "," -f1,4 | sed "1d;s/,/\t/g" | sort -k2nr > symtest.txt

  #input the symtest threshod 
  echo -e '\n'
  echo "Read the file symtest.txt in the output folder to determine the p-value threshold of symmetry tests, usually 0.01~0.1. Higher values indicate that more loci rejecting SRH hypotheses would be removed."
  echo -e '\n'
  read -p "Please input the p-value threshod:      " THRESHOLD_SYMTEST

  cat symtest.txt | grep -P "\t""0" | sort -k2nr > symtest.temp.txt
  cat symtest.txt | grep -P -v "\t""0" | cut -f1 > loci.symtest.delete
  TOTAL_LINE=$(cat symtest.temp.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p symtest.temp.txt | cut -f1)
      num=$(sed -n "$line"p symtest.temp.txt | cut -f2)
      diff=$(echo "scale=4;($num-$THRESHOLD_SYMTEST)"|bc)
      num1=`echo "$diff > 0" |bc`
      test "$num1" = 0 && echo $loci >> loci.symtest.delete
    done
  
  ls symtest | grep -v -f loci.symtest.delete > loci.symtest
  cat loci.symtest.delete | $DIR_PARALLEL/parallel -j $THREADS rm symtest/{}

  echo -e '\n'
  echo -e ""$(cat loci.symtest.delete | wc -l)" loci are removed with "$(cat loci.symtest | wc -l)" loci remaining."
  echo "The file loci.symtest contains all the file list of remaining alignments."
  echo "All the kept alignments are deposited in the folder 'symtest'."
  rm all* loci.symtest.delete iqtree symtest.temp.txt


else
  #calculate basic statistics for alignments
  mkdir -p $DIR_OUTPUT/length_site

  cp $DIR_INPUT/* $DIR_OUTPUT/length_site/
  ls $DIR_INPUT > $DIR_OUTPUT/loci.list
  cd $DIR_OUTPUT/length_site
  echo -e '\n'
  echo "Calculating alignment statistics ......"

  ln -s $DIR_PHYKIT/phykit .
  LENGTH_fun() {
    len_site=$(./phykit pis $1)
    echo -e $1"\t"$len_site > $1.len_site
  }
  export -f LENGTH_fun
  cat ../loci.list | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 1 LENGTH_fun %
  cat *.len_site | sort -k3n -k2n > ../length_site.txt
  rm *.len_site phykit

  #input the threshod 
  echo -e '\n'
  echo "Read the file length_site.txt in the output folder to determine the threshold. The 2nd~4th columns represent the number of parsimony-informative sites, the alignment length, the percentage of parsimony-informative sitesin the alignment."
  echo -e '\n'
  read -p "Please input the threshod (the minimum value):      " THRESHOLD_LEN_SITE

  TOTAL_LINE=$(cat ../length_site.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p ../length_site.txt | awk '{print $1}')
      if [ "$FILTER_METHOD" == "1" ]; then
        num=$(sed -n "$line"p ../length_site.txt | awk '{print $3}')
      elif [ "$FILTER_METHOD" == "2" ]; then
        num=$(sed -n "$line"p ../length_site.txt | awk '{print $2}')
      else
        num=$(sed -n "$line"p ../length_site.txt | awk '{print $4}')
      fi
      diff=$(echo "scale=4;($num-$THRESHOLD_LEN_SITE)"|bc)
      num1=`echo "$diff < 0" |bc`
      test "$num1" = 0 && echo $loci >> ../loci.length_site
    done

  cat ../loci.list | grep -v -f ../loci.length_site > ../loci.length_site.delete
  cat ../loci.length_site.delete | $DIR_PARALLEL/parallel -j $THREADS rm {}
  cd ..
  rm loci.list loci.length_site.delete

  echo -e '\n'
  echo "The file loci.length_site contains all the file list of remaining alignments."
  echo "All the kept alignments are deposited in the folder 'length_site'."

fi


