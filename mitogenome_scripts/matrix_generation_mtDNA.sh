#!/bin/bash
#2019.08.08 by ZF
#2023.06.13 by DSY, saturation sites detection

#copy all the MitoZ annotation results, i.e. folders *.result, into the same folder
#type 'bash matrix_generation_mtDNA.sh FOLDER_MITOZ_ANNOTATION_RESULTS'
#tools used in this script: seqkit, csvtk, parallel, FASconCAT, mafft, trimal, TransDecoder, PhyKIT, and IQ-TREE


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


#Check the threads can be used
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
until [ $THREADS -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done


#copy the busco results to 0-raw_busco/
cp -r $1 0-raw_mitoz

#Generate a file containing species list
ls 0-raw_mitoz/ | sed 's/.result//g' > species.list

#generate the gene list
cat 0-raw_mitoz/*.result/*.cds | grep ">" | cut -d ";" -f2 | sort -u > cds.list
cat 0-raw_mitoz/*.result/*.rrna | grep ">" | cut -d ";" -f2 | sort | uniq > rrna.list


#generate the raw loci sequences for each CDS/rRNA gene and and each species
#CDS
for gene in $(cat cds.list)
  do
    mkdir -p 1-raw_loci/cds-0/$gene
    for species in $(cat species.list)
      do
        cat 0-raw_mitoz/$species.result/*cds | seqkit grep -r -p "$gene" -w 0 | sed "s/R/A/g"| sed "s/Y/C/g"| sed "s/M/A/g"| sed "s/K/G/g"| sed "s/S/G/g"| sed "s/W/A/g"| sed "1c >$species" > 1-raw_loci/cds-0/$gene/$species.fna
        test -s 1-raw_loci/cds-0/$gene/$species.fna && echo "$gene sequences of $species generated..." || rm 1-raw_loci/cds-0/$gene/$species.fna
      done
  done
#rRNA
for gene in $(cat rrna.list)
  do
    mkdir -p 1-raw_loci/rrna/$gene
    for species in $(cat species.list)
      do
        cat 0-raw_mitoz/$species.result/*rrna | seqkit grep -r -p "$gene" | seqkit seq -w 0 | sed "1c >$species" > 1-raw_loci/rrna/$gene/$species.fasta
        test -s 1-raw_loci/rrna/$gene/$species.fasta && echo "$gene sequences of $species generated..." || rm 1-raw_loci/rrna/$gene/$species.fasta
      done
  done

#generate the correct coding regions of CDS
cp -r 1-raw_loci/cds-0 1-raw_loci/cds_nuc
for gene in $(cat cds.list)
  do
    cd 1-raw_loci/cds_nuc/$gene
    for fna in *.fna
      do
        TransDecoder.LongOrfs -t $fna -m 10 -G Mitochondrial-Invertebrates
        cp *dir/longest_orfs.cds $fna
        cp *dir/longest_orfs.pep $(ls $fna | cut -d . -f1).faa
        rm -rf *dir* *cmds
      done
    sed -i -n '1,2p' *
    for species in $(cat ../../../species.list); do sed -i "1c >$species" $species.*; done
    cd ../../..
  done

for gene in $(cat cds.list)
  do
    mkdir -p 1-raw_loci/cds_aa/$gene
    mv 1-raw_loci/cds_nuc/$gene/*faa 1-raw_loci/cds_aa/$gene
    sed -i "s/*//g" 1-raw_loci/cds_aa/$gene/*.faa
  done

#Merge sequences of the same locus into the fasta files
for gene in $(cat cds.list)
  do
    cat 1-raw_loci/cds_nuc/$gene/* > 1-raw_loci/cds_nuc/$gene.fna
    cat 1-raw_loci/cds_aa/$gene/* > 1-raw_loci/cds_aa/$gene.faa
  done

for gene in $(cat rrna.list)
  do
    cat 1-raw_loci/rrna/$gene/* > 1-raw_loci/rrna/$gene.fasta
  done


#Align all the faa and fna files
mkdir 2-loci_align
for gene in $(cat cds.list)
  do
    linsi --thread $THREADS 1-raw_loci/cds_aa/$gene.faa > 2-loci_align/$gene.mafft.faa
    test -s 2-loci_align/$gene.mafft.faa && echo "loci $gene has been aligned" || mafft --thread $THREADS 1-raw_loci/cds_aa/$gene.faa > 2-loci_align/$gene.mafft.faa
  done

for gene in $(cat rrna.list)
  do
    linsi --thread $THREADS 1-raw_loci/rrna/$gene.fasta > 2-loci_align/$gene.mafft.fasta
    test -s 2-loci_align/$gene.mafft.fasta && echo "loci $gene has been aligned" || mafft --thread $THREADS 1-raw_loci/rrna/$gene.fasta > 2-loci_align/$gene.mafft.fasta
  done


#Trim the alignments
mkdir -p 3-loci_trim/fna 3-loci_trim/faa 3-loci_trim/rrna
for gene in $(cat cds.list)
do
  echo "trim protein sequence of loci $gene ......"
  trimal -in 2-loci_align/$gene.mafft.faa -out 3-loci_trim/faa/$gene.aa.fas -gappyout
  echo -e '\n'
  echo "trim nucleotide sequence of loci $gene ......"
  trimal -in 2-loci_align/$gene.mafft.faa -out 3-loci_trim/fna/$gene.nuc.fas -gappyout -backtrans 1-raw_loci/cds_nuc/$gene.fna -ignorestopcodon
  echo -e '\n'
done

for gene in $(cat rrna.list)
  do
    trimal -in 2-loci_align/$gene.mafft.fasta -out 3-loci_trim/rrna/$gene.fas -gappyout
  done


#Saturate the loci
mkdir -p 4-loci_saturation/fna 4-loci_saturation/faa 4-loci_saturation/rrna 4-loci_saturation/tree
cp 3-loci_trim/fna/* 4-loci_saturation/tree/
cp 3-loci_trim/rrna/* 4-loci_saturation/tree/

cd 4-loci_saturation/tree
TREE_fun() {
    iqtree -s $1.* -m GTR+I+G -B 1000 -T 1
    }
  export -f TREE_fun
  cat ../../cds.list | parallel -I% -j 1 --max-args 1 TREE_fun %
for loci in $(cat ../../cds.list)
do
  pk_sat -a $loci.*.fas -t $loci.*.treefile >> saturation1
done
paste ../../cds.list saturation1 > tmp1.txt

TREE_fun() {
    iqtree -s $1.* -m GTR+I+G -B 1000 -T 1
    }
  export -f TREE_fun
  cat ../../rrna.list | parallel -I% -j 1 --max-args 1 TREE_fun %
for loci in $(cat ../../rrna.list)
do
  pk_sat -a $loci.fas -t $loci.fas.treefile >> saturation2
done
paste ../../rrna.list saturation2 > tmp2.txt
cat tmp1.txt tmp2.txt > saturation1.txt

cat saturation1.txt | sort -k2nr > ../saturation.txt
awk '{print $2}' ../saturation.txt > temp2
cat -b temp2 | sed "s/ //g" > temp
sed -i '1i\Loci\tValue_of_saturation' temp
csvtk -t plot line temp -x "Loci" -y "Value_of_saturation" --y-min 0 --x-max 16 --format pdf > ../saturation.pdf
cd ..

#input the value of saturation
  echo -e '\n'
  echo "Read the file saturation.txt in the output folder to determine the minimum the value of saturation. Loci of higher values of saturation are thought to be desirable."
  echo -e '\n'
  read -p "Please input the value of saturation:      " THRESHOLD_SATUR

TOTAL_LINE=$(cat saturation.txt | wc -l)
for line in $(seq $TOTAL_LINE)
  do
    loci=$(sed -n "$line"p saturation.txt | cut -f1)
    num=$(sed -n "$line"p saturation.txt | cut -f2)
    diff=$(echo "scale=4;($num-$THRESHOLD_SATUR)"|bc)
    num1=`echo "$diff < 0" |bc`
    test "$num1" = 0 && echo $loci >> loci.saturation
  done

rm -rf tree/

for loci in $(cat loci.saturation)
do
  cp ../3-loci_trim/fna/$loci.* fna/
  cp ../3-loci_trim/faa/$loci.* faa/
  cp ../3-loci_trim/rrna/$loci.* rrna/
done
cd ..


#Concatenate all the nucleotide/protein alignments in phylip format and generate partition files
mkdir -p 5-loci_concat/cds_fna 5-loci_concat/cds12_fna 5-loci_concat/cds_faa 5-loci_concat/cds_rrna 5-loci_concat/cds12_rrna 

cd 5-loci_concat/cds_faa
cp ../../4-loci_saturation/faa/* ./
perl ~/install/FASconCAT-G*/FASconCAT-G*pl -a -p -p -s -l
rm *fas

cd ../cds_fna
cp ../../4-loci_saturation/fna/* ./
perl ~/install/FASconCAT-G*/FASconCAT-G*pl -a -p -p -s -l
rm *fas

cd ../cds12_fna
cp ../../4-loci_saturation/fna/* ./
perl ~/install/FASconCAT-G*/FASconCAT-G*pl -a -p -p -s -l -d
rm *fas

cd ../cds_rrna
cp ../../4-loci_saturation/fna/* ./ && cp ../../4-loci_saturation/rrna/* ./
perl ~/install/FASconCAT-G*/FASconCAT-G*pl -a -p -p -s -l
rm *fas

cd ../cds12_rrna
cp ../cds12_fna/*phy ./cds.phy && cp ../../4-loci_saturation/rrna/* ./
perl ~/install/FASconCAT-G*/FASconCAT-G*pl -a -p -p -s -l
cp ../cds12_fna/*partition.txt partition.txt
tail -n 2 FcC_supermatrix_partition.txt >> partition.txt
rm *fas cds.phy FcC_supermatrix_partition.txt
