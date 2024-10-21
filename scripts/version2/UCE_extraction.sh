#!/bin/bash 
#2021.08.06 by ZF
#2023.07.13 by DSY
#2023.10.13 by ZF


#Type 'bash UCE_extraction.sh'
#Prepare a probe set file (e.g. Phthiraptera-2.8Kv1.fasta)
#Modify all genome assemblies endding with .fa ("SPECIES_NAME.fa") and copy them to the same folder.
#Tools faToTwoBit and parallel, as well as PHYLUCE conda enviorment, are used and will be automatically checked prior to formal analyses in this script.


##check the packages
echo "Checking the package dependency......" | tee -a log.txt
echo -e "\n" >> log.txt

#check faToTwoBit
if [ $(which faToTwoBit) ]
    then
      echo "faToTwoBit ...... OK" | tee -a log.txt
      EXE_FATOTWOBIT=$(which faToTwoBit)
      DIR_FATOTWOBIT=${EXE_FATOTWOBIT%/*}
    else
      until [ -x ${DIR_FATOTWOBIT}/faToTwoBit ]
        do
          read -p "faToTwoBit is not found. Please input its installation directory (absolute path, e.g. /usr/local/bin):      " DIR_FATOTWOBIT_TEMP
          DIR_FATOTWOBIT=$(realpath $(echo ${DIR_FATOTWOBIT_TEMP} | sed "s/'//g"))
        done
      echo "faToTwoBit ...... OK" | tee -a log.txt
fi


#check parallel
if [ $(which parallel) ]
    then
      echo "parallel ...... OK"
      EXE_PARALLEL=$(which parallel)
      DIR_PARALLEL=${EXE_PARALLEL%/*}
    else
      until [ -x ${DIR_PARALLEL}/parallel ]
        do
          read -p "parallel is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/TransDecoder-v5.5.0):      " DIR_PARALLEL_TEMP
          DIR_PARALLEL=$(realpath $(echo ${DIR_PARALLEL_TEMP} | sed "s/'//g"))
        done
      echo "PARALLEL ...... OK"
fi


#Check the threads can be used
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
until [ ${THREADS} -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done
echo "THREADS=${THREADS}" >> parameters.cfg


#input the name input directory containing all genome assemblies
read -p "Please input the name of input directory containing all genome assemblies endding with .fa ("SPECIES_NAME.fa"), e.g. /PATH/assemblies:      " DIR_ASSEMBLY_TEMP
DIR_ASSEMBLY=$(realpath $(echo ${DIR_ASSEMBLY_TEMP} | sed "s/'//g"))


#input the probe file
read -p "Please input probe file of the FASTA format, e.g. /home/zf/Desktop/materials/test/Collembola-1.8Kv1.fasta:      " INPUT_PROBE_TEMP
INPUT_PROBE=$(realpath $(echo ${INPUT_PROBE_TEMP} | sed "s/'//g"))


#copy probe set to the working folder
echo -e "\n" | tee -a log.txt
echo "copy probe set to the working folder..." | tee -a log.txt
cat ${INPUT_PROBE} > probe.fasta
echo

#copy the assemblies to 0-genomes/
echo -e "\n" | tee -a log.txt
echo "copy the assemblies to 0-genomes/..." | tee -a log.txt
cp -r ${DIR_ASSEMBLY} 0-genomes
echo


#Change the first letter of the assembly name in the folder 0-genomes/ into upper case
cd 0-genomes/
a=`ls | xargs -n1`
echo "$a" | while read line
do
  b=`echo $line | cut -c1`
  c=`echo $line |cut -c2-`
  d=`echo $b |grep -c '[a-z]'`
  if [ $d -eq 1 ]
    then
     e=`echo $b | tr 'a-z' 'A-Z'`
     mv $b$c $e$c
  fi
done
cd ..


#prepare a species list
echo "prepare a species list..." | tee -a log.txt
ls 0-genomes/ | sed "s/.fa//g" > species.list
echo


SPECIES_NAME=$(cat species.list)


#Generate a python script simplify_headers.py to simplify the sequence head
echo "simplify the sequence head..." | tee -a log.txt
echo '#!/usr/bin/python' >> simplify_headers.py
echo 'from Bio import SeqIO' >> simplify_headers.py
for SPECIES in ${SPECIES_NAME}
  do
    echo -e '\n' >> simplify_headers.py
    echo "# ${SPECIES}" >> simplify_headers.py
    echo "with open(\"${SPECIES}.fa\", "'"rU"'") as infile:" >> simplify_headers.py
    echo "  with open(\"${SPECIES}.fasta\", \"w\") as outf:" >> simplify_headers.py
    echo "    for seq in SeqIO.parse(infile, 'fasta'):" >> simplify_headers.py
    echo "      seq.name = \"\"" >> simplify_headers.py
    echo "      seq.description = \"\"" >> simplify_headers.py
    echo "      outf.write(seq.format('fasta'))" >> simplify_headers.py
  done

#Simplify the FASTA heads
cd 0-genomes/
python ../simplify_headers.py
rm *.fa
echo


#Move each genome assembly to its own directory
for SPECIES in *; do mkdir ${SPECIES%.*}; mv ${SPECIES} ${SPECIES%.*}; done

#Convert genomes to 2bit format
echo "Convert genomes to 2bit format..." | tee -a ../log.txt
cat ../species.list | ${DIR_PARALLEL}/parallel -j ${THREADS} ${DIR_FATOTWOBIT}/faToTwoBit {}/{}.fasta {}/{}.2bit
#for SPECIES in *; do $DIR_FATOTWOBIT/faToTwoBit $SPECIES/$SPECIES.fasta $SPECIES/${SPECIES%.*}.2bit; done
echo


#Generate the configure file for the genome data
echo "Generate the configure file for the genome data..." | tee -a ../log.txt
echo '[scaffolds]' >> ../target-genome.conf
for SPECIES in ${SPECIES_NAME}; do echo "${SPECIES}:../0-genomes/${SPECIES}/${SPECIES}.2bit" >> ../target-genome.conf; done


#Align bait set to the extant genome sequences
echo -e "\n" | tee -a ../log.txt
echo "Align bait set to the extant genome sequences..." | tee -a ../log.txt
cd .. && mkdir 1-uce_align && cd 1-uce_align/
phyluce_probe_run_multiple_lastzs_sqlite --db uces.sqlite --output target-genome-lastz --probefile ../probe.fasta --scaffoldlist $(cat ../species.list) --genome-base-path ../0-genomes --identity 50 --cores ${THREADS}


#Extracting FASTA sequence matching UCE loci from genome sequences
echo -e "\n" | tee -a ../log.txt
echo "Extracting FASTA sequence matching UCE loci from genome sequences..." | tee -a ../log.txt
phyluce_probe_slice_sequence_from_genomes --conf ../target-genome.conf --lastz target-genome-lastz --output target-genome-fasta --flank 400 --name-pattern "probe.fasta_v_{}.lastz.clean" | tee -a ../log.txt


#prepare the configure file for species list
echo '[all]' >> ../taxon-sets.conf
cat ../species.list >> ../taxon-sets.conf


#Match contigs to baits
echo -e "\n" | tee -a ../log.txt
phyluce_assembly_match_contigs_to_probes --contigs target-genome-fasta --probes ../probe.fasta --output in-silico-lastz --min-coverage 67 | tee -a ../log.txt


#Get match counts
echo -e "\n" | tee -a ../log.txt
cd .. && mkdir 2-uces && cd 2-uces/
phyluce_assembly_get_match_counts --locus-db ../1-uce_align/in-silico-lastz/probe.matches.sqlite --taxon-list-config ../taxon-sets.conf --taxon-group 'all' --output insilico-incomplete.conf --incomplete-matrix | tee -a ../log.txt


#Change the first letter of the species name in the folder ../../target-genome-fasta into upper case
cd ../1-uce_align/target-genome-fasta/
a=`ls | xargs -n1`
echo "$a" |while read line
do
  b=`echo $line | cut -c1`
  c=`echo $line |cut -c2-`
  d=`echo $b |grep -c '[a-z]'`
  if [ $d -eq 1 ]
    then
     e=`echo $b | tr 'a-z' 'A-Z'`
     mv $b$c $e$c
  fi
done


#Extract FASTA information of UCEs
cd ../../2-uces/
echo -e "\n" | tee -a ../log.txt
echo "Extract final UCEs..." | tee -a ../log.txt
phyluce_assembly_get_fastas_from_match_counts --contigs ../1-uce_align/target-genome-fasta --locus-db ../1-uce_align/in-silico-lastz/probe.matches.sqlite --match-count-output insilico-incomplete.conf --output insilico-incomplete.fasta --incomplete-matrix insilico-incomplete.incomplete | tee -a ../log.txt


#Generate loci list
cd .. && mkdir 3-raw_loci
sed -n '/uce-/p' 2-uces/insilico-incomplete.conf > loci.list

#Generate the sequence file for each locus
LOCI_NAME=$(cat loci.list)
cd 3-raw_loci/
for LOCI in ${LOCI_NAME}
  do
    for SPECIES in $(cat ../species.list)
      do
        echo "${LOCI}_${SPECIES} |${LOCI}" >> ${LOCI}.seq_name
      done
  done

FUN_fun() {
  awk 'FNR==NR{a[$0];next} /^>/{val=$0;sub(/^>/,"",val);flag=val in a?1:0} flag' $1.seq_name ../2-uces/insilico-incomplete.fasta > $1.fa
}
export -f FUN_fun
cat ../loci.list | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 FUN_fun %
#for LOCI in $LOCI_NAME
#  do
#    awk 'FNR==NR{a[$0];next} /^>/{val=$0;sub(/^>/,"",val);flag=val in a?1:0} flag' $LOCI.seq_name ../2-uces/insilico-incomplete.fasta > $LOCI.fa
#  done
rm *name


#Simplify the sequence heads
for LOCI in ${LOCI_NAME}; do sed -i "s/${LOCI}_\| |${LOCI}//g" ${LOCI}.fa; done


#Summarize the loci number for each species
echo -e "\n" | tee -a ../log.txt
SPECIES_NAME=$(cat ../species.list)
for SPECIES in ${SPECIES_NAME}
do
  for LOCI in ${LOCI_NAME}
  do
  if [ -z $(grep -o ${SPECIES} ${LOCI}.fa) ]
    then
      echo "$LOCI in $SPECIES does not exist" >> log.txt
  fi
  done
done


#Filter loci having too few taxa (less than three)
cd .. && mkdir 4-loci_filter
TOTAL_TAXA=$(cat species.list | wc -l)
for LOCI in ${LOCI_NAME}
do
  if [ $(grep -o $LOCI 3-raw_loci/log.txt | wc -l) -lt `expr ${TOTAL_TAXA} - 2` ]
    then
        cp 3-raw_loci/${LOCI}.fa 4-loci_filter
        echo "${LOCI}" >> list.uce
  fi 
done
