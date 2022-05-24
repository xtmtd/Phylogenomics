#Estimate divergence time using MCMCTree for large AA dataset
#Type "sh mcmctree_AA.sh IPQTREE_PARTITION_FILE FOLDER_alignments", e.g., sh mcmctree_AA.sh FcC_supermatrix_partition.txt.best_model.nex  4-loci_trim/faa_bmge
#'IPQTREE_PARTITION_FILE' is the model merging strategy estimated from IQ-TREE




#check FASconCAT
until [ -s $DIR_FASconCAT/FASconCAT*.pl ]
    do
      read -p "DIR_FASconCAT is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/FASconCAT-G-1.04):      " DIR_FASconCAT
    done
echo "FASconCAT ...... OK"



#check readAL
if [ $(which readal) ]
    then
      echo "readal ...... OK"
      EXE_READAL=$(which readal)
      DIR_READAL=${EXE_READAL%/*}
    else
      until [ -x $DIR_READAL/readal ]
        do
          read -p "readal is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/trimal-1.4.1/source):      " DIR_READAL
        done
      echo "readal ...... OK"
fi


#check mcmctree
if [ $(which mcmctree) ]
    then
      echo "mcmctree ...... OK"
      EXE_MCMCTREE=$(which mcmctree)
      DIR_MCMCTREE=${EXE_MCMCTREE%/*}
    else
      until [ -x $DIR_MCMCTREE/mcmctree ]
        do
          read -p "mcmctree is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/paml-4.9j/bin):      " DIR_MCMCTREE
        done
      echo "mcmctree ...... OK"
fi


#check codeml
if [ $(which codeml) ]
    then
      echo "codeml ...... OK"
      EXE_CODEML=$(which codeml)
      DIR_CODEML=${EXE_CODEML%/*}
    else
      until [ -x $DIR_CODEML/codeml ]
        do
          read -p "codeml is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/paml-4.9j/bin):      " DIR_CODEML
        done
      echo "codeml ...... OK"
fi


#check the number of samples for mcmctree
read -p "Please input the number of samples/generations kept for MCMCTree (e.g. 10000, 20000):      " NSAMPLE
  until [ $NSAMPLE -gt 0 ]
    do
      read -p "Please input the correct integer for the number of samples/generations kept for MCMCTree (e.g. 10000, 20000):      " NSAMPLE
    done


#check parallel
if [ $(which parallel) ]
    then
      echo "parallel ...... OK"
      EXE_PARALLEL=$(which parallel)
      DIR_PARALLEL=${EXE_PARALLEL%/*}
    else
      until [ -x $DIR_PARALLEL/parallel ]
        do
          read -p "parallel is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/TransDecoder-v5.5.0):      " DIR_PARALLEL
        done
      echo "PARALLEL ...... OK"
fi


#Check the threads can be used
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
  until [ $THREADS -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done



mkdir -p 0-partition 1-phylip 2-BV 3-divergence/run1 3-divergence/run2

#generate the phylip sequence file compatible with PAML
cat $1 | grep "charset" | cut -d "=" -f 1 | cut -d " " -f 4 > temp
for id in $(seq "$(cat temp | wc -l)")
  do
    mkdir par && cd par
    cat ../temp | sed -n ""$id"p" | tr "_" "\n" | cut -d "." -f 1 > ../par.list
    for loci in $(cat ../par.list); do cp $2/"$loci"* ./; done
    perl $DIR_FASconCAT/FASconCAT-G*pl -s -l
    mv FcC*fas ../0-partition/p$id.fas
    cd ..
    rm -rf par*
  done

for id in $(seq "$(cat temp | wc -l)"); do readal -in 0-partition/p$id.fas -out 1-phylip/p$id.phy -phylip_paml; done

cat 1-phylip/*phy > input.phy


#modify mcmctree.ctl file
DIR_CUR=$(pwd)
cd $DIR_MCMCTREE && cd ..
cp ./mcmctree.ctl $DIR_CUR/
cd $DIR_CUR

sed -i "2c seqfile = $PWD/input.phy" mcmctree.ctl
sed -i "3c treefile = $PWD/input.tre" mcmctree.ctl
NUM_LOCI=$(ls 0-partition/ | wc -l)
sed -i "s/ndata = 3/ndata = "$NUM_LOCI"/g" mcmctree.ctl
sed -i "s/seqtype = 0/seqtype = 2/g" mcmctree.ctl
sed -i "s/usedata = 1/usedata = 3/g" mcmctree.ctl
sed -i "s/clock = 3/clock = 2/g" mcmctree.ctl
sed -i "s/RootAge = <1.0/RootAge = /g" mcmctree.ctl
sed -i "s/alpha = 0/alpha = 0.5/g" mcmctree.ctl
sed -i "s/cleandata = 0/cleandata = 1/g" mcmctree.ctl
sed -i "s/BDparas = 1 1 0/BDparas = 1 1 0.1/g" mcmctree.ctl
sed -i "s/rgene_gamma = 2 2/rgene_gamma = 2 20 1/g" mcmctree.ctl
sed -i "s/sigma2_gamma = 1 10/sigma2_gamma = 1 10 1/g" mcmctree.ctl
sed -i "s/finetune = 1: 0.1  0.1  0.1  0.01 .5/finetune = 0: .1  .1  .1  .1 .1 .1/g" mcmctree.ctl
sed -i "s/burnin = 2000/burnin = 10000/g" mcmctree.ctl
sed -i "s/sampfreq = 2/sampfreq = 5/g" mcmctree.ctl
sed -i "s/nsample = 20000/nsample = "$NSAMPLE"/g" mcmctree.ctl


#generate in.BV
cd 2-BV
cp ../mcmctree.ctl ./
mcmctree mcmctree.ctl

rm out.BV rst*
cd $DIR_MCMCTREE && cd ..
cp ./dat/lg.dat $DIR_CUR/2-BV/
cd $DIR_CUR/2-BV

for file in tmp*ctl; do sed -i "s/model = 0/model = 2/g" $file; sed -i "s/aaRatefile = /aaRatefile = lg.dat/g" $file; sed -i "s/method = 0/method = 1/g" $file; done
#for file in tmp*ctl; do $DIR_CODEML/codeml $file; cat rst2 >> in.BV; done
export DIR_CODEML
CODEML_fun() {
	$DIR_CODEML/codeml $1
	cat rst2 >> in.BV
}
export -f CODEML_fun
find . -name "tmp*ctl" | parallel -j $THREADS CODEML_fun

#rm l* out r* SeedUsed t*
cd ..


#estimate divergence time
for run in 1 2;
  do
    cp 2-BV/in.BV 3-divergence/run$run/
    cp mcmctree.ctl 3-divergence/run$run/
    sed -i "s/usedata = 3/usedata = 2/g" 3-divergence/run$run/mcmctree.ctl
  done

cd 3-divergence/run1
nohup mcmctree mcmctree.ctl &

cd ../run2
nohup mcmctree mcmctree.ctl &








