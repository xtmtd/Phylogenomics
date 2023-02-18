#!/bin/bash
#2021.08.04 wrote by ZF
#2023.01.13 revised by DSY

#filter loci by detecting alignment length, number of parsimony-informative sites, percentage of parsimony-informative sites in the alignment, GC content, compositional heterogeneity (RCV, Relative composition variability), evolutionary rate (average pairwise identity) and likelihood mapping criteria.
#Type 'bash loci_filtering_alignment-based.sh'
#parallel, PhyKIT, seqkit, csvtk and IQ-TREE may be required


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


#check Seqkit
if [ $(which seqkit) ]
    then
      echo "SEQKIT ...... OK"
      EXE_SEQKIT=$(which seqkit)
      DIR_SEQKIT=${EXE_SEQKIT%/*}
    else
      until [ -x $DIR_SEQKIT/seqkit ]
  do
    read -p "SEQKIT is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_SEQKIT_TEMP
    DIR_SEQKIT=$(realpath $(echo $DIR_SEQKIT_TEMP | sed "s/'//g"))
  done
      echo "SEQKIT ...... OK"
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

  awk '{print $2}' ../rcv.txt > ../temp2
  cat -b ../temp2 | sed "s/ //g" > ../temp
  sed -i '1i Loci_number\tRCV' ../temp
  csvtk -t plot line ../temp -x Loci_number -y RCV --format pdf > ../RCV.pdf
  rm ../temp ../temp2

  #input the RCV threshold 
  echo -e '\n'
  echo "Read the file rcv.txt in the output folder to determine the RCV threshold. Lower RCV values are thought to be desirable because they represent a lower composition bias in an alignment."
  echo -e '\n'
  read -p "Please input the RCV threshold:      " THRESHOLD_RCV

  TOTAL_LINE=$(cat ../rcv.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p ../rcv.txt | cut -f1)
      num=$(sed -n "$line"p ../rcv.txt | cut -f2)
      diff=$(echo "scale=4;($num-$THRESHOLD_RCV)"|bc)
      num1=`echo "$diff >= 0" |bc`
      test "$num1" = 0 && echo $loci >> ../loci.rcv
    done
  
  cat ../loci.rcv | sed 's/[ ]//g' > ../loci.rcv1
  cat ../loci.rcv1 | sed 's/[\t]//g' > ../loci.rcv2
  cat ../loci.list | sed 's/[ ]//g' > ../loci.list1
  cat ../loci.list1 | sed 's/[\t]//g' > ../loci.list2
  rm -rf ../loci.rcv ../loci.list ../loci.rcv1 ../loci.list1
  mv ../loci.rcv2 ../loci.rcv
  mv ../loci.list2 ../loci.list
  cp ../loci.list ../loci.rcv.delete
  for id in $(cat ../loci.rcv)
  do
    sed -i "/^$id/d" ../loci.rcv.delete
  done
  cat ../loci.rcv.delete | $DIR_PARALLEL/parallel -j $THREADS rm {}
  cd ..
  rm loci.list loci.rcv.delete

  cd $DIR_OUTPUT/rcv/
  $DIR_PHYKIT/phykit cc -a ../loci.rcv -p rcv
  rm *.occupancy
  mv rcv* ..
  AVE=$(echo "scale=2;$(cat $DIR_OUTPUT/rcv.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2)/$(cat $DIR_OUTPUT/loci.rcv | wc -l)"|bc)

  echo -e '\n'
  echo "RCV with $THRESHOLD_RCV (threshold) has preserved $(cat $DIR_OUTPUT/rcv.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/loci.rcv | wc -l) loci (the average site of each locus is $AVE)." | tee -a $DIR_OUTPUT/summary.statistics
  rm ../rcv.partition ../rcv.fa
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

  awk '{print $2}' ../GC.txt > ../temp2
  cat -b ../temp2 | sed "s/ //g" > ../temp
  sed -i '1i\Loci_number\tGC' ../temp
  csvtk -t plot line ../temp -x Loci_number -y GC --format pdf > ../GC.pdf
  rm ../temp ../temp2

  #input the GC threshold 
  echo -e '\n'
  echo "Read the file GC.txt in the output folder to determine the GC percentage threshold. Very high (often > 0.8) GC values are thought to be unreliable."
  echo -e '\n'
  read -p "Please input the GC percentage threshold:      " THRESHOLD_GC

  TOTAL_LINE=$(cat ../GC.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p ../GC.txt | cut -f1)
      num=$(sed -n "$line"p ../GC.txt | cut -f2)
      diff=$(echo "scale=4;($num-$THRESHOLD_GC)"|bc)
      num1=`echo "$diff >= 0" |bc`
      test "$num1" = 0 && echo $loci >> ../loci.GC
    done

  cat ../loci.GC | sed 's/[ ]//g' > ../loci.GC1
  cat ../loci.GC1 | sed 's/[\t]//g' > ../loci.GC2
  cat ../loci.list | sed 's/[ ]//g' > ../loci.list1
  cat ../loci.list1 | sed 's/[\t]//g' > ../loci.list2
  rm -rf ../loci.GC ../loci.list ../loci.GC1 ../loci.list1
  mv ../loci.GC2 ../loci.GC
  mv ../loci.list2 ../loci.list
  cp ../loci.list ../loci.GC.delete
  for id in $(cat ../loci.GC)
  do
    sed -i "/^$id/d" ../loci.GC.delete
  done
  cat ../loci.GC.delete | $DIR_PARALLEL/parallel -j $THREADS rm {}
  cd ..
  rm loci.list loci.GC.delete

  cd $DIR_OUTPUT/GC/
  $DIR_PHYKIT/phykit cc -a ../loci.GC -p GC
  rm *.occupancy
  mv GC* ..
  AVE=$(echo "scale=2;$(cat $DIR_OUTPUT/GC.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2)/$(cat $DIR_OUTPUT/loci.GC | wc -l)"|bc)

  echo -e '\n'
  echo "GC with $THRESHOLD_GC (threshold) has preserved $(cat $DIR_OUTPUT/GC.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/loci.GC | wc -l) loci (the average site of each locus is $AVE)." | tee -a $DIR_OUTPUT/summary.statistics
  rm ../GC.partition ../GC.fa
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

  awk '{print $2}' ../evolutionary_rate.txt > ../temp2
  cat -b ../temp2 | sed "s/ //g" > ../temp
  sed -i '1i\Loci_number\tEvolutionary_rate' ../temp
  csvtk -t plot line ../temp -x Loci_number -y Evolutionary_rate --format pdf > ../Evolutionary_rate.pdf
  rm ../temp ../temp2

  #input the evolutionary rate threshold 
  echo -e '\n'
  echo "Read the file evolutionary_rate.txt in the output folder to determine the rate threshold. Very high values indicate that those 'fast' genes may be more suitable for shallow phylogenies."
  echo -e '\n'
  read -p "Please input the evolutionary_rate threshold:      " THRESHOLD_RATE

  TOTAL_LINE=$(cat ../evolutionary_rate.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p ../evolutionary_rate.txt | cut -f1)
      num=$(sed -n "$line"p ../evolutionary_rate.txt | cut -f2)
      diff=$(echo "scale=4;($num-$THRESHOLD_RATE)"|bc)
      num1=`echo "$diff >= 0" |bc`
      test "$num1" = 0 && echo $loci >> ../loci.evolutionary_rate
    done

  cat ../loci.evolutionary_rate | sed 's/[ ]//g' > ../loci.evolutionary_rate1
  cat ../loci.evolutionary_rate1 | sed 's/[\t]//g' > ../loci.evolutionary_rate2
  cat ../loci.list | sed 's/[ ]//g' > ../loci.list1
  cat ../loci.list1 | sed 's/[\t]//g' > ../loci.list2
  rm -rf ../loci.evolutionary_rate ../loci.list ../loci.evolutionary_rate1 ../loci.list1
  mv ../loci.evolutionary_rate2 ../loci.evolutionary_rate
  mv ../loci.list2 ../loci.list
  cp ../loci.list ../loci.evolutionary_rate.delete
  for id in $(cat ../loci.evolutionary_rate)
  do
    sed -i "/^$id/d" ../loci.evolutionary_rate.delete
  done
  cat ../loci.evolutionary_rate.delete | $DIR_PARALLEL/parallel -j $THREADS rm {}
  cd ..
  rm loci.list loci.evolutionary_rate.delete

  cd $DIR_OUTPUT/evolutionary_rate
  $DIR_PHYKIT/phykit cc -a ../loci.evolutionary_rate -p evolutionary_rate
  rm *.occupancy
  mv evolutionary_rate* ..
  AVE=$(echo "scale=2;$(cat $DIR_OUTPUT/evolutionary_rate.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2)/$(cat $DIR_OUTPUT/loci.evolutionary_rate | wc -l)"|bc)

  echo -e '\n'
  echo "Evolutionary rate with $THRESHOLD_RATE (threshold) has preserved $(cat $DIR_OUTPUT/evolutionary_rate.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/loci.evolutionary_rate | wc -l) loci (the average site of each locus is $AVE)." | tee -a $DIR_OUTPUT/summary.statistics
  rm ../evolutionary_rate.partition ../evolutionary_rate.fa
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
  cat *.quartet | sort -k2nr > ../likelihood_mapping_quartet.txt
  mv *.iqtree ../iqtree/
  cd .. 
  rm -rf tempLM

  awk '{print $2}' likelihood_mapping_quartet.txt > temp2
  cat -b temp2 | sed "s/ //g" > temp
  sed -i '1i\Loci_number\tLikelihood_mapping_quartet' temp
  csvtk -t plot line temp -x Loci_number -y Likelihood_mapping_quartet --format pdf > Likelihood_mapping.pdf
  rm temp temp2

  #input the quartet threshold 
  echo -e '\n'
  echo "Read the file likelihood_mapping_quartet.txt in the output folder to determine the quartet threshold. High values indicate that strong phylogenetic informativeness/signal."
  echo -e '\n'
  read -p "Please input the quartet threshold:      " THRESHOLD_QUARTET

  TOTAL_LINE=$(cat likelihood_mapping_quartet.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p likelihood_mapping_quartet.txt | cut -f1)
      num=$(sed -n "$line"p likelihood_mapping_quartet.txt | cut -f2)
      diff=$(echo "scale=2;($num-$THRESHOLD_QUARTET)"|bc)
      num1=`echo "$diff >= 0" |bc`
      test "$num1" = 1 && echo $loci >> loci.likelihood_mapping
    done

  cat loci.likelihood_mapping | $DIR_PARALLEL/parallel -j $THREADS cp $DIR_INPUT/{} likelihood_mapping/

  cd $DIR_OUTPUT/likelihood_mapping
  $DIR_PHYKIT/phykit cc -a ../loci.likelihood_mapping -p likelihood_mapping
  rm *.occupancy
  mv likelihood_mapping* ..
  AVE=$(echo "scale=2;$(cat $DIR_OUTPUT/likelihood_mapping.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2)/$(cat $DIR_OUTPUT/loci.likelihood_mapping | wc -l)"|bc)

  echo -e '\n'
  echo "Likelihood mapping with $THRESHOLD_QUARTET (quartet threshold) has preserved $(cat $DIR_OUTPUT/likelihood_mapping.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/loci.likelihood_mapping | wc -l) loci (the average site of each locus is $AVE)." | tee -a $DIR_OUTPUT/summary.statistics
  rm ../likelihood_mapping.partition ../likelihood_mapping.fa
  echo "The file loci.likelihood_mapping contains all the file list of remaining alignments."
  echo "All the kept alignments are deposited in the folder 'likelihood_mapping.'"


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

  awk '{print $2}' symtest.txt > temp2
  cat -b temp2 | sed "s/ //g" > temp
  sed -i '1i\Loci_number\tSymtest' temp
  csvtk -t plot line temp -x Loci_number -y Symtest --format pdf > Symtest.pdf
  rm temp temp2

  #input the symtest threshold 
  echo -e '\n'
  echo "Read the file symtest.txt in the output folder to determine the p-value threshold of symmetry tests, usually 0.01~0.1. Higher values indicate that more loci rejecting SRH hypotheses would be removed."
  echo -e '\n'
  read -p "Please input the p-value threshold:      " THRESHOLD_SYMTEST

  cat symtest.txt | grep -P "\t""0" | sort -k2nr > symtest.temp.txt
  cat symtest.txt | grep -P -v "\t""0" | cut -f1 > loci.symtest.delete
  TOTAL_LINE=$(cat symtest.temp.txt | wc -l)
  for line in $(seq $TOTAL_LINE)
    do
      loci=$(sed -n "$line"p symtest.temp.txt | cut -f1)
      num=$(sed -n "$line"p symtest.temp.txt | cut -f2)
      diff=$(echo "scale=4;($num-$THRESHOLD_SYMTEST)"| bc)
      num1=`echo "$diff > 0" | bc`
      test "$num1" = 0 && echo $loci >> loci.symtest.delete
    done

  ls symtest/ > loci
  cat loci | sed 's/[ ]//g' > loci1
  cat loci1 | sed 's/[\t]//g' > loci2
  cat loci.symtest.delete | sed 's/[ ]//g' > loci.symtest.delete1
  cat loci.symtest.delete1 | sed 's/[\t]//g' > loci.symtest.delete2
  rm -rf loci.symtest.delete1 loci.symtest.delete loci loci1
  mv loci2 loci.symtest
  mv loci.symtest.delete2 loci.symtest.delete
  for id in $(cat loci.symtest.delete)
  do
    sed -i "/^$id/d" loci.symtest
  done
  cat loci.symtest.delete | $DIR_PARALLEL/parallel -j $THREADS rm symtest/{}

  cd $DIR_OUTPUT/symtest
  $DIR_PHYKIT/phykit cc -a ../loci.symtest -p symtest
  rm *.occupancy
  mv symtest* ..
  AVE=$(echo "scale=2;$(cat $DIR_OUTPUT/symtest.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2)/$(cat $DIR_OUTPUT/loci.symtest | wc -l)"|bc)
  cd ..

  echo -e '\n'
  echo -e ""$(cat loci.symtest.delete | wc -l)" loci are removed, and "$(cat loci.symtest | wc -l)" loci are remained with $(cat $DIR_OUTPUT/symtest.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites (the average site of each locus is $AVE)." | tee -a $DIR_OUTPUT/summary.statistics
  echo "The file loci.symtest contains all the file list of remaining alignments."
  echo "All the kept alignments are deposited in the folder 'symtest'."
  rm all* loci.symtest.delete iqtree symtest.temp.txt symtest.fa symtest.partition


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

  if [ "$FILTER_METHOD" == "2" ]; then
  awk '{print $2}' ../length_site.txt > ../temp1
  cat ../temp1 | csvtk -t -H sort -k1:n > ../temp2
  cat -b ../temp2 | sed "s/ //g" > ../temp
  sed -i '1i\Loci_number\tParsimony_sites' ../temp
  csvtk -t plot line ../temp -x Loci_number -y Parsimony_sites --format pdf > ../Parsimony_sites.pdf
  rm ../temp ../temp2 ../temp1
  elif [ "$FILTER_METHOD" == "1" ]; then
  awk '{print $3}' ../length_site.txt > ../temp1
  cat ../temp1 | csvtk -t -H sort -k1:n > ../temp2
  cat -b ../temp2 | sed "s/ //g" > ../temp
  sed -i '1i\Loci_number\tAlignment_length' ../temp
  csvtk -t plot line ../temp -x Loci_number -y Alignment_length --format pdf > ../Alignment_length.pdf
  rm ../temp ../temp2 ../temp1
  else
  awk '{print $4}' ../length_site.txt > ../temp1
  cat ../temp1 | csvtk -t -H sort -k1:n > ../temp2
  cat -b ../temp2 | sed "s/ //g" > ../temp
  sed -i '1i\Loci_number\tPercentage_of_parsimony_sites' ../temp
  csvtk -t plot line ../temp -x Loci_number -y Percentage_of_parsimony_sites --format pdf > ../Percentage_of_parsimony_sites.pdf
  rm ../temp ../temp2 ../temp1
  fi

  #input the threshold 
  echo -e '\n'
  echo "Read the file length_site.txt in the output folder to determine the threshold. The 2nd~4th columns represent the number of parsimony-informative sites, the alignment length, the percentage of parsimony-informative sites in the alignment."
  echo -e '\n'
  read -p "Please input the threshold (the minimum value):      " THRESHOLD_LEN_SITE

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

  cat ../loci.length_site | sed 's/[ ]//g' > ../loci.length_site1
  cat ../loci.length_site1 | sed 's/[\t]//g' > ../loci.length_site2
  cat ../loci.list | sed 's/[ ]//g' > ../loci.list1
  cat ../loci.list1 | sed 's/[\t]//g' > ../loci.list2
  rm -rf ../loci.length_site ../loci.list ../loci.length_site1 ../loci.list1
  mv ../loci.length_site2 ../loci.length_site
  mv ../loci.list2 ../loci.list
  cp ../loci.list ../loci.length_site.delete
  for id in $(cat ../loci.length_site)
  do
    sed -i "/^$id/d" ../loci.length_site.delete
  done
  cat ../loci.length_site.delete | $DIR_PARALLEL/parallel -j $THREADS rm {}
  cd ..
  rm loci.list loci.length_site.delete

  cd $DIR_OUTPUT/length_site
  $DIR_PHYKIT/phykit cc -a ../loci.length_site -p length_site
  rm *.occupancy
  mv length_site* ..
  AVE=$(echo "scale=2;$(cat $DIR_OUTPUT/length_site.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2)/$(cat $DIR_OUTPUT/loci.length_site | wc -l)"|bc)

  echo -e '\n'
  echo "Length/parsimony-informative sites with $THRESHOLD_LEN_SITE (threshold) has preserved $(cat $DIR_OUTPUT/length_site.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/loci.length_site | wc -l) loci (the average site of each locus is $AVE)." | tee -a $DIR_OUTPUT/summary.statistics
  rm ../length_site.partition ../length_site.fa
  echo "The file loci.length_site contains all the file list of remaining alignments."
  echo "All the kept alignments are deposited in the folder 'length_site'."

fi


