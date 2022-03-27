#generate the supermaxtrix, partition and occupancy for alignments
#Type 'sh matrix_generation.sh'
#parallel and PhyKIT is required


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
          DIR_PHYKIT=$(echo $DIR_PHYKIT_TEMP | sed "s/'//g")
        done
      echo "PhyKIT ...... OK"
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


#input the prefix for the generated matrix-related files
read -p "Please input the name of a prefix for the generated matrix-related files, e.g. DATASET1:      " PREFIX


#input the taxon occupancy values for the matrices
read -p "Please input the minimum percentage value for taxa occupancy, usually ranging from 50% to 100%, e.g. 50, 75, 90:      " OCCUPANCY


cd $DIR_INPUT
ls > ../$PREFIX.alignments
phykit cc -a ../$PREFIX.alignments -p $PREFIX
rm ../$PREFIX.alignments

mkdir -p $DIR_OUTPUT/all $DIR_OUTPUT/matrix$OCCUPANCY/alignments
mv $PREFIX.* $DIR_OUTPUT/all

cd $DIR_OUTPUT/all 
TOTAL_LINE=$(cat $PREFIX.occupancy | wc -l)
for line in $(seq $TOTAL_LINE)
  do
    loci=$(sed -n "$line"p $PREFIX.occupancy | awk '{print $1}')
    num=$(sed -n "$line"p $PREFIX.occupancy | awk '{print $4}')
    OCCU=$(echo "scale=4;($OCCUPANCY/100)"|bc)
    diff=$(echo "scale=4;($num-$OCCU)"|bc)
    num1=`echo "$diff < 0" |bc`
    test "$num1" = 0 && echo $loci >> $DIR_OUTPUT/matrix$OCCUPANCY/$PREFIX.alignments
  done

cd $DIR_OUTPUT/matrix$OCCUPANCY/alignments
for file in $(cat ../$PREFIX.alignments); do cp $DIR_INPUT/$file ./; done
phykit cc -a ../$PREFIX.alignments -p $PREFIX
mv  $PREFIX.* $DIR_OUTPUT/matrix$OCCUPANCY/
rm -rf $DIR_OUTPUT/all

echo "Individual loci alignments, concatenated matrix and partition file are deposited in the OUTPUT/matrix$OCCUPANCY."
echo "$OCCUPANCY% occupancy matrix has $(cat $DIR_OUTPUT/matrix$OCCUPANCY/$PREFIX.partition | tail -n -1 | cut -d "=" -f2 | cut -d "-" -f2) sites and $(cat $DIR_OUTPUT/matrix$OCCUPANCY/$PREFIX.partition | wc -l) loci." | tee -a $DIR_OUTPUT/summary.matrices
cd $DIR_CURR
