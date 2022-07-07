#!/bin/bash
#2021.10.10 wrote by DSY
#2022.07.06 revised by DSY

#detection of outlier long branches in collections of phylogenetic trees
#type "bash treeshrink.sh"
#TreeShrink, python and R is required


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

#input the name of input directory containing all gene trees
read -p "Please input the name of input directory containing all gene trees, e.g. gene_trees:      " DIR_INPUT_TREE_TEMP
DIR_INPUT_TREE_TEMP1=$(echo $DIR_INPUT_TREE_TEMP | sed "s/'//g")
cd $DIR_INPUT_TREE_TEMP1 && DIR_INPUT_TREE=$(echo $PWD) && cd $DIR_CURR

#input the name of input directory containing all loci alignments
read -p "Please input the name of input directory containing all loci alignments, e.g. matrix75/alignment:      " DIR_INPUT_ALIGN_TEMP
DIR_INPUT_ALIGN_TEMP1=$(echo $DIR_INPUT_ALIGN_TEMP | sed "s/'//g")
cd $DIR_INPUT_ALIGN_TEMP1 && DIR_INPUT_ALIGN=$(echo $PWD) && cd $DIR_CURR

#input the name of output directory
read -p "Please input the name of output directory, or an existing directory:      " DIR_OUTPUT_TEMP
DIR_OUTPUT_TEMP1=$(echo $DIR_OUTPUT_TEMP | sed "s/'//g")
test -d $DIR_OUTPUT_TEMP1 || mkdir -p $DIR_OUTPUT_TEMP1
cd $DIR_OUTPUT_TEMP1 && DIR_OUTPUT=$(echo $PWD) && cd $DIR_CURR

#input the number of α threshold (α = 0.05 is default)
read -p "Please input the number of α threshold. Default is 0.05 (e.g. 0.05):      " THRESHOLD
  until [ $THRESHOLD ]
    do
      read -p "Please input the number of α threshold. Default is 0.05 (e.g. 0.05):      " THRESHOLD
    done


#Check the number of job and the threads per job
read -p "Please input the number of threads/cores used for each TreeShrink analysis (e.g. 8):      " THREADS
  until [ $THREADS -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores used for each TreeShrink analysis (e.g. 8):      " THREADS
    done


#copy the alignments and gene trees, and extract the loci
ls $DIR_INPUT_ALIGN/ > $DIR_OUTPUT/loci.list
mkdir -p $DIR_OUTPUT/treeshrink && cd $DIR_OUTPUT/treeshrink
cat ../loci.list | $DIR_PARALLEL/parallel -j $THREADS mkdir {}
cat ../loci.list | $DIR_PARALLEL/parallel -j $THREADS cp $DIR_INPUT_ALIGN/{} {}
cat ../loci.list | $DIR_PARALLEL/parallel -j $THREADS cp $DIR_INPUT_TREE/{}.treefile {}
cd ..

#treeshrink analysis
TREESHRINK_fun() {
  run_treeshrink.py -i treeshrink/ -t $1.treefile -q "$THRESHOLD" -a $1 -O $1.treeshrink
}
export -f TREESHRINK_fun
cat loci.list | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 1 TREESHRINK_fun %

#extract filtered sequences
for loci in $(cat loci.list)
do
  [[ `cat treeshrink/$loci/*.txt | wc -l` -eq 0 ]] && echo "$loci" >> list.keep
done

  mv loci.list list.clean
  for id in $(cat list.keep)
  do
    sed -i "/^$id/d" list.clean
  done

mkdir treeshrink_keep
cat list.keep | $DIR_PARALLEL/parallel -j $THREADS mv treeshrink/{}/{} treeshrink_keep/
rm -rf treeshrink/

cat list.clean | wc -l > temp
if [ "temp" == "0" ]
then
mkdir treeshrink_clean summary_txt 
cat list.clean | $DIR_PARALLEL/parallel -j $THREADS mv treeshrink/{}/{}.treeshrink.fas treeshrink_fas/
cat list.clean | $DIR_PARALLEL/parallel -j $THREADS mv treeshrink/{}/*.txt summary_txt/
rm -rf temp
echo -e '\n'
echo "Filtered alignments (removed the taxa with abnormally long branches) have been deposited in the folder $DIR_OUTPUT/treeshrink_clean/. All these alignments should be redone the gene trees."
echo "List of loci containing outlier sequences was saved in the file $DIR_OUTPUT/list.clean"
echo "Alignments without abnormally long branches have been deposited in the folder $DIR_OUTPUT/treeshrink_keep/."
echo "The taxa of each alignment must more than 10. If not, these alignments will not be analysed. Check them carefully!"
else
rm -rf list.clean temp
echo -e '\n'
echo "Congratulations! All alignments have no abnormally long branches."
echo -e '\n'
echo "Alignments without abnormally long branches have been deposited in the folder $DIR_OUTPUT/treeshrink_keep/."
echo "The taxa of each alignment must more than 10. If not, these alignments will not be analysed. Check them carefully!"
fi
