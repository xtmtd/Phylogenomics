#generate individual gene trees for a set of alignments
#Type 'sh gene_trees.sh'
#IQ-TREE v2 and parallel are required


##Checking the package dependency
echo "Checking the package dependency......"

#check IQ-TREE
if [ $(which iqtree) ]
    then
      echo "IQ-TREE ...... OK"
      EXE_IQTREE=$(which iqtree)
      DIR_IQTREE=${EXE_IQTREE%/*}
    else
      until [ -x $EXE_IQTREE ]
        do
          read -p "IQ-TREE is not found. Please input its location (absolute path, e.g. /usr/bin/iqtree):      " EXE_IQTREE_TEMP
          EXE_IQTREE=$(echo $EXE_IQTREE_TEMP | sed "s/'//g")
        done
      echo "IQ-TREE ...... OK"
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
          DIR_PARALLEL=$(echo $DIR_PARALLEL_TEMP | sed "s/'//g")
        done
      echo "PARALLEL ...... OK"
fi


DIR_CURR=$(echo $PWD)

#input the name of input directory
read -p "Please input the name of input directory containing all alignments, e.g. 4-trim/clipkit-kpi:      " DIR_INPUT_TEMP
DIR_INPUT_TEMP1=$(echo $DIR_INPUT_TEMP | sed "s/'//g")
cd $DIR_INPUT_TEMP1 && DIR_INPUT=$(echo $PWD) && cd $DIR_CURR


#input the name of output directory
read -p "Please input the name of output directory, or an existing directory:      " DIR_OUTPUT_TEMP
DIR_OUTPUT_TEMP1=$(echo $DIR_OUTPUT_TEMP | sed "s/'//g")
test -d $DIR_OUTPUT_TEMP1 || mkdir -p $DIR_OUTPUT_TEMP1
cd $DIR_OUTPUT_TEMP1 && DIR_OUTPUT=$(echo $PWD) && cd $DIR_CURR



#Check the type of input alignments
read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide   " ALIGN_TYPE
  until [ $ALIGN_TYPE -gt 0 -a $ALIGN_TYPE -lt 3 ]
    do
      read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide      " ALIGN_TYPE
    done


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

else
  read -p "Please input the option for DNA substitution model: 1. automatically ModelFinder searching; 2. restricted to HKY/GTR model to reduce computational burden (-mset HKY,GTR)         " MODEL_INPUT
    until [ $MODEL_INPUT -gt 0 -a $MODEL_INPUT -lt 3 ]
      do
        read -p "Please input the option for DNA substitution model: 1. automatically ModelFinder searching; 2. restricted to HKY/GTR model to reduce computational burden (-mset HKY,GTR)         " MODEL_INPUT
      done
  [ "$MODEL_INPUT" == "1" ] && MODEL=""
  [ "$MODEL_INPUT" == "2" ] && MODEL="--mset HKY,GTR"

fi


#Check the number of job and the threads per job
echo -e "\n"
echo "Tasks of tree calculations will be parallelized: product of the number of IQ-TREE jobs (-j) and the number of used threads/cores (-T) per job cannot be greater than the number of physical threads/cores"

read -p "Please input the number of threads/cores used for each IQ-TREE analysis (e.g. 8):      " THREADS
  until [ $THREADS -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores used for each IQ-TREE analysis (e.g. 8):      " THREADS
    done

read -p "Please input the number of IQ-TREE jobs/tasks (e.g. 4):      " JOBS
  until [ $THREADS -gt 0 ]
    do
      read -p "Please input the number of IQ-TREE jobs/tasks (e.g. 4)::      " JOBS
    done

echo -e "\n"
echo "Calculating gene trees ......"


cp -r $DIR_INPUT TEMP_trees && cd TEMP_trees
mkdir -p $DIR_OUTPUT/gene_trees
ls $DIR_INPUT > ../temp.list

cat ../temp.list | parallel -j $JOBS $EXE_IQTREE -s {} $MODEL -B 1000 -T $THREADS
mv *treefile $DIR_OUTPUT/gene_trees/
cd ..
rm -rf TEMP_trees temp.list

echo -e "\n"
echo "All genes trees have been deposited in OUTPUT_DIRECTORY/gene_trees."





