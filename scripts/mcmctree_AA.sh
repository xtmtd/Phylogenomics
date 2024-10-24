#!/bin/bash
#2021.08.04 by ZF
#2022.08.03 by ZF, optimization of in.BV generation

#Estimate divergence time using MCMCTree for large AA dataset with approximate likelihood computation applied
#Type "bash mcmctree_AA.sh"
#A tree with fossil calibrations, a partition file estimated from IQ-TREE (such as XXXbest_model.nex), and a folder containing all alignments are required
#Phykit, PAML, trimal, parallel and csvtk are required for this script



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


#check readAL
if [ $(which readal) ]
    then
      echo "readal ...... OK"
      EXE_READAL=$(which readal)
      DIR_READAL=${EXE_READAL%/*}
    else
      until [ -x $DIR_READAL/readal ]
        do
          read -p "readal (a module of trimAl package) is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/trimal-1.4.1/source):      " DIR_READAL_TEMP
          DIR_READAL=$(realpath $(echo $DIR_READAL_TEMP | sed "s/'//g"))
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
          read -p "mcmctree (PAML module) is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/paml-4.9j/bin):      " DIR_MCMCTREE_TEMP
          DIR_MCMCTREE=$(realpath $(echo $DIR_MCMCTREE_TEMP | sed "s/'//g"))
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
          read -p "codeml (PAML module) is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/paml-4.9j/bin):      " DIR_CODEML_TEMP
          DIR_CODEML=$(realpath $(echo $DIR_CODEML_TEMP | sed "s/'//g"))
        done
      echo "codeml ...... OK"
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


#check the number of samples for mcmctree
read -p "Default burnin and sample frequency (sampfreq) are set as 50000 and 5 in this script. Please input the number of samples/generations kept for MCMCTree, 10000~20000 could be enough for most cases):      " NSAMPLE
  until [ $NSAMPLE -gt 0 ]
    do
      read -p "Default burnin and sample frequency (sampfreq) are set as 50000 and 5 in this script. Please input the number of samples/generations kept for MCMCTree, 10000~20000 could be enough for most cases):      " NSAMPLE
    done


#input the name input directory containing all alignments
read -p "Please input the name of input directory containing all alignments, e.g. /PATH/alignments:      " DIR_ALIGNMENT_TEMP
DIR_ALIGNMENT=$(realpath $(echo $DIR_ALIGNMENT_TEMP | sed "s/'//g"))


#input the partition file generated by iqtree 
read -p "Please input the partition file 'XXXbest_model.nex' (with its path) generated by IQ-TREE, e.g. /home/zf/Desktop/materials/test/test.partition.best_model.nex:      " INPUT_PARTITION_TEMP
INPUT_PARTITION=$(realpath $(echo $INPUT_PARTITION_TEMP | sed "s/'//g"))


#input the tree file containing fossil calibrations 
read -p "Please input the tree file containing fossil calibrations, e.g. /home/zf/Desktop/materials/test/input.tre. The first line conatins the number of taxa and the number of trees, respectively, e.g. '6  1', following by newick rooted tree format with calibrations (unit 100 Mya) appended on nodes, e.g (Zootermopsis_nevadensis,((Pediculus_humanus,(Thrips_palmi,Maconellicoccus_hirsutus)'>3.146<3.827'),(Hypothenemus_hampei,Drosophila_melanogaster)'>2.955<3.589'))'<4.192';. Root age must be constrained in the tree:      " INPUT_TREE_TEMP
INPUT_TREE=$(realpath $(echo $INPUT_TREE_TEMP | sed "s/'//g"))


DIR_CURR=$(echo $PWD)
mkdir -p 0-partition 1-phylip 2-BV 3-divergence/run1 3-divergence/run2

#generate the phylip sequence file compatible with PAML
cat $INPUT_PARTITION | grep "charset" | cut -d "=" -f 1 | cut -d " " -f 4 > temp
for id in $(seq -w "$(cat temp | wc -l)")
  do
    mkdir par_mcmctree && cd par_mcmctree
    cat ../temp | sed -n ""$id"p" | tr "_" "\n" > ../par.id
    for loci in $(cat ../par.id); do cp $DIR_ALIGNMENT/"$loci"* .; done
    ls > ../par.alignments
    $DIR_PHYKIT/phykit cc -a ../par.alignments -p p$id
    mv p*fa ../0-partition/p$id.fas
    cd ..
    rm -rf par.id par.alignments par_mcmctree
  done


#transform fasta files into paml_phylip format
for id in $(seq -w "$(cat temp | wc -l)"); do $DIR_READAL/readal -in 0-partition/p$id.fas -out 1-phylip/p$id.phy -phylip_paml; done
cat 1-phylip/*phy > input.phy

rm temp


#modify mcmctree.ctl file
cd $DIR_MCMCTREE && cd ..
cp ./mcmctree.ctl $DIR_CURR/
cd $DIR_CURR

sed -i "2c seqfile = $DIR_CURR/input.phy" mcmctree.ctl
sed -i "3c treefile = $INPUT_TREE" mcmctree.ctl
NUM_LOCI=$(ls 0-partition/ | wc -l)
sed -i "s/ndata = 3/ndata = "$NUM_LOCI"/g" mcmctree.ctl
sed -i "s/seqtype = 0/seqtype = 2/g" mcmctree.ctl
sed -i "s/usedata = 1/usedata = 3/g" mcmctree.ctl
sed -i "s/clock = 3/clock = 2/g" mcmctree.ctl
sed -i "s/RootAge = <1.0/RootAge = /g" mcmctree.ctl
sed -i "s/cleandata = 0/cleandata = 1/g" mcmctree.ctl
sed -i "s/BDparas = 1 1 0/BDparas = 1 1 0.1/g" mcmctree.ctl
sed -i "s/rgene_gamma = 2 2/rgene_gamma = 2 20 1/g" mcmctree.ctl
sed -i "s/sigma2_gamma = 1 10/sigma2_gamma = 1 10 1/g" mcmctree.ctl
sed -i "s/finetune = 1: 0.1  0.1  0.1  0.01 .5/finetune = 0: .1  .1  .1  .1 .1 .1/g" mcmctree.ctl
sed -i "s/burnin = 2000/burnin = 50000/g" mcmctree.ctl
sed -i "s/sampfreq = 2/sampfreq = 5/g" mcmctree.ctl
sed -i "s/nsample = 20000/nsample = "$NSAMPLE"/g" mcmctree.ctl


#generate in.BV
cd 2-BV
ls ../1-phylip/*.phy | cut -d "/" -f3 | cut -d "." -f1 > list.partition

cd $DIR_MCMCTREE && cd ../dat/
cp lg.dat $DIR_CURR/2-BV/
cd $DIR_CURR/2-BV
ln -s $DIR_CODEML/codeml .
ln -s $DIR_MCMCTREE/mcmctree .

BV_fun() {
  mkdir $1 && cd $1
  cp ../lg.dat .
  cp ../../mcmctree.ctl .
  sed -i "2c seqfile = ../../1-phylip/$1.phy" mcmctree.ctl
  sed -i "6c ndata = 1" mcmctree.ctl
  timeout 10s ../mcmctree mcmctree.ctl
  sed -i "s/model = 0/model = 2/g" tmp0001.ctl
  sed -i "s/aaRatefile = /aaRatefile = lg.dat/g" tmp0001.ctl
  ../codeml tmp0001.ctl
  cd ..  
}
export -f BV_fun
cat list.partition | $DIR_PARALLEL/parallel -I% -j $THREADS --max-args 1 BV_fun %

for folder in $(cat list.partition); do cat $folder/rst2 >> in.BV; done
rm -rf p* l* codeml mcmctree


#estimate divergence time
cd $DIR_CURR/3-divergence
for id in 1 2
  do
    cp $DIR_CURR/2-BV/in.BV run$id
    cp $DIR_CURR/mcmctree.ctl  run$id
    sed -i "s/usedata = 3/usedata = 2/g" run*/mcmctree.ctl
    ln -s $DIR_MCMCTREE/mcmctree run$id/mcmctree
  done

MCMCTREE_RUN_fun() {
  cd $1
  ./mcmctree mcmctree.ctl > log.txt
  rm mcmctree
  cd ..
}
export -f MCMCTREE_RUN_fun
find . -name "run*" | $DIR_PARALLEL/parallel -I% -j 2 --max-args 1 MCMCTREE_RUN_fun %
echo

echo "You can view the posterior tree 'FigTree.tre' (in the folders runXX) using tool FigTree."
echo


#generate convergence plot to check the convergence between two runs
#obtain posterior times
cd $DIR_CURR/3-divergence

cat run1/out | grep "^t_" | awk '{print $1}' > t1
cat run1/out | grep "^t_" | awk '{print $2,$4,$5}' | sed "s/,//g;s/)//g;s/ /\t/g" > t2
cat run1/out | grep "^t_" | awk '{print $2,$4,$5}' | sed "s/,//g;s/)//g;s/ /\t/g" | awk '{print ($3-$2)}' > t3
cat run2/out | grep "^t_" | awk '{print $2,$4,$5}' | sed "s/,//g;s/)//g;s/ /\t/g" > t4
cat run2/out | grep "^t_" | awk '{print $2,$4,$5}' | sed "s/,//g;s/)//g;s/ /\t/g" | awk '{print ($3-$2)}' > t5
paste t* | sed "1i node\tmean_posterior_run1\tlower_95%CI_run1\tupper_95%CI_run1\tCI_width_run1\tmean_posterior_run2\tlower_95%CI_run2\tupper_95%CI_run2\tCI_width_run2" > posterior_times.txt
rm t*


#compare the posterior means of each node between two runs
#calculate pearson correlation between two runs
cat posterior_times.txt | $DIR_CSVTK/csvtk -t corr -T -f mean_posterior_run1,mean_posterior_run2 1>>Pearson_correlation_coefficient.txt 2>&1
CORR=$(cat Pearson_correlation_coefficient.txt | cut -f3)
$DIR_CSVTK/csvtk -t plot line posterior_times.txt -x mean_posterior_run1 -y mean_posterior_run2 --format pdf -o convergence_plot_posterior.pdf

echo -e "Pearson correlation coefficient of the posterior means between two runs is "$CORR". You can view the convergence plot 'convergence_plot_posterior.pdf' in the folder 3-divergence to check the convergence. The points should fall almost perfectly on the straight line; if not, increase the generation number of MCMC runs."
echo

#generate infinite-sites plot assess whether collecting additional molecular data would improve the analysis
cat posterior_times.txt | $DIR_CSVTK/csvtk -t corr -T -f mean_posterior_run1,CI_width_run1 1>>Pearson_correlation_coefficient.txt 2>&1
$DIR_CSVTK/csvtk -t plot line posterior_times.txt -x mean_posterior_run1 -y CI_width_run1 --format pdf -o infinite-sites_plot_run1_posterior.pdf
cat posterior_times.txt | $DIR_CSVTK/csvtk -t corr -T -f mean_posterior_run2,CI_width_run2 1>>Pearson_correlation_coefficient.txt 2>&1
$DIR_CSVTK/csvtk -t plot line posterior_times.txt -x mean_posterior_run2 -y CI_width_run2 --format pdf -o infinite-sites_plot_run2_posterior.pdf

echo "View the infinite-sites plots of each run for posterior distribution. It helps to assess whether collecting additional molecular data would improve the analysis. As the number of sites and loci tend to infinity, a plot of mean posterior times vs. credibility interval widths will tend to a straight line. Pearson correlation coefficient of two variables are also calculated and kept in the file 'Pearson_correlation_coefficient.txt'"
echo

#compare prior times to posterior times
#change "usedata = 0" and run mcmctree again to obtain priot times
mkdir -p run1/prior_distribution run2/prior_distribution
for id in 1 2
  do
    cd run$id/prior_distribution
    mv ../in.BV .
    cp ../mcmctree.ctl .
    SEED=$(cat ../SeedUsed)
    sed -i "s/seed = -1/seed = "$SEED"/g" mcmctree.ctl
    sed -i "s/usedata = 2/usedata = 0/g" mcmctree.ctl
    ln -s $DIR_MCMCTREE/mcmctree mcmctree
    cd ../..
  done

PRIOR_RUN_fun() {
  cd $1/prior_distribution
  mcmctree mcmctree.ctl > log.txt
  rm mcmctree in.BV
  cd ../..
}
export -f PRIOR_RUN_fun
find . -name "run*" | $DIR_PARALLEL/parallel -I% -j 2 --max-args 1 PRIOR_RUN_fun %

echo "See the priot time estimates in the folders runXX/prior_divergence."
echo

#obtain prior times
cd $DIR_CURR/3-divergence

cat run1/prior_distribution/out | grep "^t_" | awk '{print $1}' > t1
cat run1/prior_distribution/out | grep "^t_" | awk '{print $2,$4,$5}' | sed "s/,//g;s/)//g;s/ /\t/g" > t2
cat run1/prior_distribution/out | grep "^t_" | awk '{print $2,$4,$5}' | sed "s/,//g;s/)//g;s/ /\t/g" | awk '{print ($3-$2)}' > t3
cat run2/prior_distribution/out | grep "^t_" | awk '{print $2,$4,$5}' | sed "s/,//g;s/)//g;s/ /\t/g" > t4
cat run2/prior_distribution/out | grep "^t_" | awk '{print $2,$4,$5}' | sed "s/,//g;s/)//g;s/ /\t/g" | awk '{print ($3-$2)}' > t5
paste t* | sed "1i node\tmean_prior_run1\tlower_95%CI_run1\tupper_95%CI_run1\tCI_width_run1\tmean_prior_run2\tlower_95%CI_run2\tupper_95%CI_run2\tCI_width_run2" > prior_times.txt
rm t*

#compare the prior means of each node between two runs
#calculate pearson correlation between two runs
cat prior_times.txt | $DIR_CSVTK/csvtk -t corr -T -f mean_prior_run1,mean_prior_run2 1>>Pearson_correlation_coefficient.txt 2>&1
CORR=$(cat Pearson_correlation_coefficient.txt | tail -n 1 | cut -f3)
$DIR_CSVTK/csvtk -t plot line prior_times.txt -x mean_prior_run1 -y mean_prior_run2 --format pdf -o convergence_plot_prior.pdf

echo -e "Pearson correlation coefficient of the prior means between two runs is "$CORR". You can view the convergence plot 'convergence_plot_prior.pdf' in the folder 3-divergence to check the convergence. The points should fall almost perfectly on the straight line; if not, increase the generation number of MCMC runs."
echo

#generate infinite-sites plot assess whether collecting additional molecular data would improve the analysis
cat prior_times.txt | $DIR_CSVTK/csvtk -t corr -T -f mean_prior_run1,CI_width_run1 1>>Pearson_correlation_coefficient.txt 2>&1
$DIR_CSVTK/csvtk -t plot line prior_times.txt -x mean_prior_run1 -y CI_width_run1 --format pdf -o infinite-sites_plot_run1_prior.pdf
cat prior_times.txt | $DIR_CSVTK/csvtk -t corr -T -f mean_prior_run2,CI_width_run2 1>>Pearson_correlation_coefficient.txt 2>&1
$DIR_CSVTK/csvtk -t plot line prior_times.txt -x mean_prior_run2 -y CI_width_run2 --format pdf -o infinite-sites_plot_run2_prior.pdf

echo "View the infinite-sites plots of each run for prior distribution. It helps to assess whether collecting additional molecular data would improve the analysis. As the number of sites and loci tend to infinity, a plot of mean prior times vs. credibility interval widths will tend to a straight line. Pearson correlation coefficient of two variables are also calculated and kept in the file 'Pearson_correlation_coefficient.txt'"
echo

#mean posterior vs mean prior
cat posterior_times.txt | cut -f2,6 > t1
cat prior_times.txt | cut -f2,6 > t2
paste t1 t2 > t3
cat t3 | $DIR_CSVTK/csvtk -t corr -T -f mean_posterior_run1,mean_prior_run1 1>>Pearson_correlation_coefficient.txt 2>&1
$DIR_CSVTK/csvtk -t plot line t3 -x mean_posterior_run1 -y mean_prior_run1 --format pdf -o posterior_vs_prior_plot_run1.pdf
cat t3 | $DIR_CSVTK/csvtk -t corr -T -f mean_posterior_run2,mean_prior_run2 1>>Pearson_correlation_coefficient.txt 2>&1
$DIR_CSVTK/csvtk -t plot line t3 -x mean_posterior_run2 -y mean_prior_run2 --format pdf -o posterior_vs_prior_plot_run2.pdf
rm t*

echo "View the mean posterior vs mean prior plots of each run. It helps to assess that whether prior distribution estimated from fossil calibrations are consistent with those posterior dustribution. Pearson correlation coefficient of two variables are also calculated and kept in the file 'Pearson_correlation_coefficient.txt'"

