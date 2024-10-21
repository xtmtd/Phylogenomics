#!/bin/bash
#2021.08.04 by ZF
#2022.08.03 by ZF, optimization of in.BV generation
#2023.05.31 by DSY, suitable for pmal-4.10, and add the nucleotide analysis
#2023.10.24 by ZF

#Estimate divergence time using MCMCTree for large AA dataset with approximate likelihood computation applied
#Type "bash mcmctree_AA.sh"
#A tree with fossil calibrations, a partition file estimated from IQ-TREE (such as XXXbest_model.nex), and a folder containing all alignments are required
#Phykit, PAML >= v4.10, trimal, parallel and csvtk are required for this script



#check Phykit
if [ $(which phykit) ]
    then
      echo "PhyKIT ...... OK"
      EXE_PHYKIT=$(which phykit)
      DIR_PHYKIT=${EXE_PHYKIT%/*}
    else
      until [ -x ${DIR_PHYKIT}/phykit ]
        do
          read -p "PhyKIT is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_PHYKIT_TEMP
          DIR_PHYKIT=$(realpath $(echo ${DIR_PHYKIT_TEMP} | sed "s/'//g"))
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
      until [ -x ${DIR_READAL}/readal ]
        do
          read -p "readal (a module of trimAl package) is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/trimal-1.4.1/source):      " DIR_READAL_TEMP
          DIR_READAL=$(realpath $(echo ${DIR_READAL_TEMP} | sed "s/'//g"))
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
      until [ -x ${DIR_MCMCTREE}/mcmctree ]
        do
          read -p "mcmctree (PAML module) is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/paml-4.10.6/bin):      " DIR_MCMCTREE_TEMP
          DIR_MCMCTREE=$(realpath $(echo ${DIR_MCMCTREE_TEMP} | sed "s/'//g"))
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
      until [ -x ${DIR_CODEML}/codeml ]
        do
          read -p "codeml (PAML module) is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/paml-4.10.6/bin):      " DIR_CODEML_TEMP
          DIR_CODEML=$(realpath $(echo ${DIR_CODEML_TEMP} | sed "s/'//g"))
        done
      echo "codeml ...... OK"
fi

#check baseml
if [ $(which baseml) ]
    then
      echo "baseml ...... OK"
      EXE_BASEML=$(which baseml)
      DIR_BASEML=${EXE_BASEML%/*}
    else
      until [ -x ${DIR_BASEML}/baseml ]
        do
          read -p "baseml (PAML module) is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/paml-4.10.6/bin):      " DIR_BASEML_TEMP
          DIR_BASEML=$(realpath $(echo ${DIR_BASEML_TEMP} | sed "s/'//g"))
        done
      echo "baseml ...... OK"
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


#check csvtk
if [ $(which csvtk) ]
    then
      echo "csvtk ...... OK"
      EXE_CSVTK=$(which csvtk)
      DIR_CSVTK=${EXE_CSVTK%/*}
    else
      until [ -x ${DIR_CSVTK}/csvtk ]
        do
          read -p "CSVTK is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_CSVTK_TEMP
          DIR_CSVTK=$(realpath $(echo ${DIR_CSVTK_TEMP} | sed "s/'//g"))
        done
      echo "csvtk ...... OK"
fi


#Check the threads can be used
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
  until [ ${THREADS} -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done

#Check the type of input alignments
read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide   " ALIGN_TYPE
  until [[ ${ALIGN_TYPE} -gt 0 && ${ALIGN_TYPE} -lt 3 ]]
    do
      read -p "Please input the option for input alignments: 1. amino acid; 2. nucleotide      " ALIGN_TYPE
    done


#check the number of samples for mcmctree
read -p "Default burnin and sample frequency (sampfreq) are set as 50000 and 5 in this script. Please input the number of samples/generations kept for MCMCTree, 10000~20000 could be enough for most cases):      " NSAMPLE
  until [ ${NSAMPLE} -gt 0 ]
    do
      read -p "Default burnin and sample frequency (sampfreq) are set as 50000 and 5 in this script. Please input the number of samples/generations kept for MCMCTree, 10000~20000 could be enough for most cases):      " NSAMPLE
    done


#input the name input directory containing all alignments
read -p "Please input the name of input directory containing all alignments, e.g. /PATH/alignments:      " DIR_ALIGNMENT_TEMP
DIR_ALIGNMENT=$(realpath $(echo ${DIR_ALIGNMENT_TEMP} | sed "s/'//g"))
until [ -d ${DIR_ALIGNMENT} ]
   do
    read -p "Please input the name of input directory containing all alignments, e.g. /PATH/alignments:      " DIR_ALIGNMENT_TEMP
    DIR_ALIGNMENT=$(realpath $(echo ${DIR_ALIGNMENT_TEMP} | sed "s/'//g"))
   done
echo "input directory containing all alignments ...... OK"


#whether to partition
read -p "Whether to partition the sequence based on the average pairwise identity (API): 1. no; 2. yes   " PAR_TYPE
  until [[ ${PAR_TYPE} -gt 0 && ${PAR_TYPE} -lt 3 ]]
    do
      read -p "Whether to partition the sequence based on the average pairwise identity (API): 1. no; 2. yes      " PAR_TYPE
    done


#input the tree file containing fossil calibrations 
read -p "Please input the tree file containing fossil calibrations, e.g. /home/zf/Desktop/materials/test/input.tre. The first line conatins the number of taxa and the number of trees, respectively, e.g. '6  1', following by newick rooted tree format with calibrations (unit 100 Mya) appended on nodes, e.g (Zootermopsis_nevadensis,((Pediculus_humanus,(Thrips_palmi,Maconellicoccus_hirsutus)'>3.146<3.827'),(Hypothenemus_hampei,Drosophila_melanogaster)'>2.955<3.589'))'<4.192';. Root age must be constrained in the tree:      " INPUT_TREE_TEMP
INPUT_TREE=$(realpath $(echo ${INPUT_TREE_TEMP} | sed "s/'//g"))
until [ -s ${INPUT_TREE} ]
   do
      read -p "Please input the tree file containing fossil calibrations, e.g. /home/zf/Desktop/materials/test/input.tre. The first line conatins the number of taxa and the number of trees, respectively, e.g. '6  1', following by newick rooted tree format with calibrations (unit 100 Mya) appended on nodes, e.g (Zootermopsis_nevadensis,((Pediculus_humanus,(Thrips_palmi,Maconellicoccus_hirsutus)'>3.146<3.827'),(Hypothenemus_hampei,Drosophila_melanogaster)'>2.955<3.589'))'<4.192';. Root age must be constrained in the tree:      " INPUT_TREE_TEMP
      INPUT_TREE=$(realpath $(echo ${INPUT_TREE_TEMP} | sed "s/'//g"))
   done
echo "input the tree file containing fossil calibrations ...... OK"


#input the name of output directory
read -p "Please input the name of output directory, or an existing directory:      " DIR_OUTPUT_TEMP
DIR_OUTPUT_TEMP1=$(echo ${DIR_OUTPUT_TEMP} | sed "s/'//g")
test -d ${DIR_OUTPUT_TEMP1} || mkdir -p ${DIR_OUTPUT_TEMP1}
cd ${DIR_OUTPUT_TEMP1} && DIR_OUTPUT=$(echo ${PWD}) && cd ${DIR_OUTPUT}


mkdir -p 0-partition 1-phylip 2-BV 3-divergence/run1 3-divergence/run2 tmp


if [ ${PAR_TYPE} -eq 2 ]; then
  ls ${DIR_ALIGNMENT} > ${DIR_OUTPUT}/loci.list   #generate loci list
  #calculate average pairwise identity
  echo "Calculating average pairwise identity ......"
  
  export DIR_ALIGNMENT DIR_PHYKIT DIR_OUTPUT
  tmp_fun() {
    echo -e "$1""\t""$(${DIR_PHYKIT}/phykit pi ${DIR_ALIGNMENT}/$1 | awk 'NR == 1 {print $2}')" > ./tmp/$1
  }
  export -f tmp_fun
  cat loci.list | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 tmp_fun %
  
  cat tmp/* | sort -k2n > API.txt
  rm -rf tmp loci.list
  
  i=$(cat API.txt | wc -l)
  j=`expr $i / 4`
  split -l $j API.txt -d -a 1 api_
  cat api_0 | awk '{print $1}' | cut -d "." -f1 | tr "\n" "_" | sed -e 's/,$/\n/' | awk '{print $0}' > api0
  cat api_1 | awk '{print $1}' | cut -d "." -f1 | tr "\n" "_" | sed -e 's/,$/\n/' | awk '{print $0}' > api1
  cat api_2 | awk '{print $1}' | cut -d "." -f1 | tr "\n" "_" | sed -e 's/,$/\n/' | awk '{print $0}' > api2
  cat api_3 | awk '{print $1}' | cut -d "." -f1 | tr "\n" "_" | sed -e 's/,$/\n/' | awk '{print $0}' > api3
  cat api0 api1 api2 api3 > temp
  rm api* API.txt

else
  ls ${DIR_ALIGNMENT} > ${DIR_OUTPUT}/loci.list   #generate loci list
  export DIR_ALIGNMENT DIR_PHYKIT DIR_OUTPUT
  cat loci.list | cut -d "." -f1 | tr "\n" "_" | sed -e 's/,$/\n/' | awk '{print $0}' > temp
  rm -rf loci.list tmp

fi

  
#construct alignments
for id in $(seq -w "$(cat temp | wc -l)")
  do
    mkdir par_mcmctree && cd par_mcmctree
    sed -n "${id}p" ../temp | tr "_" "\n" > ../par.id
    for loci in $(cat ../par.id); do cp ${DIR_ALIGNMENT}/${loci}* .; done
    ls > ../par.alignments
    ${DIR_PHYKIT}/phykit cc -a ../par.alignments -p p${id}
    mv p*fa ../0-partition/p${id}.fas
    cd ..
    rm -rf par.id par.alignments par_mcmctree
  done


#transform fasta files into paml_phylip format
for id in $(seq -w "$(cat temp | wc -l)"); do ${DIR_READAL}/readal -in 0-partition/p${id}.fas -out 1-phylip/p${id}.phy -phylip_paml; done
cat 1-phylip/*phy > input.phy
rm temp*


#modify mcmctree.ctl file
cd ${DIR_MCMCTREE} && cd ..
cp ./examples/mcmctree.ctl ${DIR_OUTPUT}/
cd ${DIR_OUTPUT}

sed -i 's/\r$//' mcmctree.ctl #transform dos format to Linux format
sed -i "/seqfile/c seqfile = ${DIR_OUTPUT}/input.phy" mcmctree.ctl
sed -i "/treefile/c treefile = ${INPUT_TREE}" mcmctree.ctl
NUM_LOCI=$(ls 0-partition/ | wc -l)
sed -i "/ndata = 3/c ndata = ${NUM_LOCI}" mcmctree.ctl
sed -i "/usedata =/c usedata = 3" mcmctree.ctl
sed -i "/clock =/c clock = 2" mcmctree.ctl
sed -i "/RootAge =/c RootAge = " mcmctree.ctl
sed -i "/cleandata =/c cleandata = 1" mcmctree.ctl
sed -i "/BDparas =/c BDparas = 1 1 0.1" mcmctree.ctl
sed -i "/rgene_gamma =/c rgene_gamma = 2 20 1" mcmctree.ctl
sed -i "/sigma2_gamma =/c sigma2_gamma = 1 10 1" mcmctree.ctl
sed -i "/finetune =/c finetune = 0: .1  .1  .1  .1 .1 .1" mcmctree.ctl
sed -i "/burnin =/c burnin = 50000" mcmctree.ctl
sed -i "/sampfreq =/c sampfreq = 5" mcmctree.ctl
sed -i "/nsample =/c nsample = ${NSAMPLE}" mcmctree.ctl


#protein mode
if [ ${ALIGN_TYPE} -eq 1 ]; then

  sed -i "/seqtype =/c seqtype = 2" mcmctree.ctl

  #generate in.BV
  cd ${DIR_OUTPUT}/1-phylip
  ls *.phy | cut -d "/" -f3 | cut -d "." -f1 > ${DIR_OUTPUT}/2-BV/list.partition

  cd ${DIR_MCMCTREE} && cd ../dat/
  cp lg.dat ${DIR_OUTPUT}/2-BV/ #using LG model for protein alignments
  cd ${DIR_OUTPUT}/2-BV

  export DIR_CODEML DIR_MCMCTREE DIR_OUTPUT
  BV_fun() {
    mkdir $1 && cd $1
    cp ${DIR_OUTPUT}/2-BV/lg.dat .
    cp ${DIR_OUTPUT}/mcmctree.ctl .
    sed -i "/seqfile =/c seqfile = ${DIR_OUTPUT}/1-phylip/"$1".phy" mcmctree.ctl
    sed -i "/^ndata =/c ndata = 1" mcmctree.ctl
    timeout 10s ${DIR_MCMCTREE}/mcmctree mcmctree.ctl
    sed -i "/model =/c model = 2" tmp0001.ctl
    sed -i "/aaRatefile =/c aaRatefile = lg.dat" tmp0001.ctl
    ${DIR_CODEML}/codeml tmp0001.ctl
    cd ..  
  }
  export -f BV_fun
  cat list.partition | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 BV_fun %

  for folder in $(cat list.partition); do cat ${folder}/rst2 >> in.BV; done
  rm -rf p* l*

#nucleotide mode
else

  sed -i "/seqtype =/c seqtype = 0" mcmctree.ctl
  sed -i "/model =/c model = 5" mcmctree.ctl #HKY+G5 model
  sed -i "/alpha =/c alpha = 0.5" mcmctree.ctl

  #generate in.BV
  cd ${DIR_OUTPUT}/1-phylip
  ls *.phy | cut -d "/" -f3 | cut -d "." -f1 > ${DIR_OUTPUT}/2-BV/list.partition
  cd ${DIR_OUTPUT}/2-BV

  export DIR_BASEML DIR_MCMCTREE DIR_OUTPUT
  BV_fun() {
    mkdir $1 && cd $1
    cp ${DIR_OUTPUT}/mcmctree.ctl .
    sed -i "/seqfile =/c seqfile = ${DIR_OUTPUT}/1-phylip/"$1".phy" mcmctree.ctl
    sed -i "/^ndata =/c ndata = 1" mcmctree.ctl
    timeout 10s ${DIR_MCMCTREE}/mcmctree mcmctree.ctl
    sed -i "/model =/c model = 4" tmp0001.ctl #HKY model
    ${DIR_BASEML}/baseml tmp0001.ctl
    cd ..  
  }
  export -f BV_fun
  cat list.partition | ${DIR_PARALLEL}/parallel -I% -j ${THREADS} --max-args 1 BV_fun %

  for folder in $(cat list.partition); do cat ${folder}/rst2 >> in.BV; done
  rm -rf p* l*

fi

#estimate divergence time
cd ${DIR_OUTPUT}/3-divergence
for id in 1 2
  do
    cp ${DIR_OUTPUT}/2-BV/in.BV run${id}/
    cp ${DIR_OUTPUT}/mcmctree.ctl  run${id}/
    sed -i "/usedata =/c usedata = 2 *" run*/mcmctree.ctl #remember to add "*" at the end, otherwise a bug will happen
  done

export DIR_MCMCTREE
MCMCTREE_RUN_fun() {
  cd $1
  ${DIR_MCMCTREE}/mcmctree mcmctree.ctl > log.txt
  cd ..
}
export -f MCMCTREE_RUN_fun
find . -name "run*" | ${DIR_PARALLEL}/parallel -I% -j 2 --max-args 1 MCMCTREE_RUN_fun %
echo

echo "You can view the posterior tree 'FigTree.tre' (in the folders runXX) using tool FigTree."
echo


#generate convergence plot to check the convergence between two runs
#obtain posterior times
cd ${DIR_OUTPUT}/3-divergence

awk '/^t_/ {print $1,$2,$4,$5}' run1/out | sed "s/,//g;s/)//g;s/ /\t/g" | awk -v OFS="\t" '{print $1, $2, $3, $4, $4-$3}' > t1
awk '/^t_/ {print $2,$4,$5}' run2/out | sed "s/,//g;s/)//g;s/ /\t/g" | awk -v OFS="\t" '{print $1, $2, $3, $3-$2}' > t2
paste t* | sed "1i node\tmean_posterior_run1\tlower_95%CI_run1\tupper_95%CI_run1\tCI_width_run1\tmean_posterior_run2\tlower_95%CI_run2\tupper_95%CI_run2\tCI_width_run2" > posterior_times.txt
rm t*


#compare the posterior means of each node between two runs
#calculate pearson correlation between two runs
cat posterior_times.txt | ${DIR_CSVTK}/csvtk -t corr -T -f mean_posterior_run1,mean_posterior_run2 1>>Pearson_correlation_coefficient.txt 2>&1
CORR=$(cut -f3 Pearson_correlation_coefficient.txt)
${DIR_CSVTK}/csvtk -t plot line posterior_times.txt -x mean_posterior_run1 -y mean_posterior_run2 --format pdf -o convergence_plot_posterior.pdf

echo -e "Pearson correlation coefficient of the posterior means between two runs is "$CORR". You can view the convergence plot 'convergence_plot_posterior.pdf' in the folder 3-divergence to check the convergence. The points should fall almost perfectly on the straight line; if not, increase the generation number of MCMC runs."
echo

#generate infinite-sites plot assess whether collecting additional molecular data would improve the analysis
cat posterior_times.txt | ${DIR_CSVTK}/csvtk -t corr -T -f mean_posterior_run1,CI_width_run1 1>>Pearson_correlation_coefficient.txt 2>&1
${DIR_CSVTK}/csvtk -t plot line posterior_times.txt -x mean_posterior_run1 -y CI_width_run1 --format pdf -o infinite-sites_plot_run1_posterior.pdf
cat posterior_times.txt | ${DIR_CSVTK}/csvtk -t corr -T -f mean_posterior_run2,CI_width_run2 1>>Pearson_correlation_coefficient.txt 2>&1
${DIR_CSVTK}/csvtk -t plot line posterior_times.txt -x mean_posterior_run2 -y CI_width_run2 --format pdf -o infinite-sites_plot_run2_posterior.pdf

echo "View the infinite-sites plots of each run for posterior distribution. It helps to assess whether collecting additional molecular data would improve the analysis. As the number of sites and loci tend to infinity, a plot of mean posterior times vs. credibility interval widths will tend to a straight line. Pearson correlation coefficient of two variables are also calculated and kept in the file 'Pearson_correlation_coefficient.txt'"
echo

#compare prior times to posterior times
#change "usedata = 0" and run mcmctree again to obtain priot times
cd ${DIR_OUTPUT}/3-divergence
mkdir -p run1/prior_distribution run2/prior_distribution
for id in 1 2
  do
    cd run${id}/prior_distribution
    mv ../in.BV .
    cp ../mcmctree.ctl .
    SEED=$(cat ../SeedUsed)
    sed -i "s/seed = -1/seed = ${SEED}/g" mcmctree.ctl
    sed -i "s/usedata = 2/usedata = 0/g" mcmctree.ctl
    cd ../..
  done

export DIR_MCMCTREE
PRIOR_RUN_fun() {
  cd $1/prior_distribution
  ${DIR_MCMCTREE}/mcmctree mcmctree.ctl > log.txt
  rm in.BV
  cd ../..
}
export -f PRIOR_RUN_fun
find . -name "run*" | ${DIR_PARALLEL}/parallel -I% -j 2 --max-args 1 PRIOR_RUN_fun %

echo "See the priot time estimates in the folders runXX/prior_divergence."
echo

#obtain prior times
cd ${DIR_OUTPUT}/3-divergence

awk '/^t_/ {print $1,$2,$4,$5}' run1/prior_distribution/out | sed "s/,//g;s/)//g;s/ /\t/g" | awk -v OFS="\t" '{print $1, $2, $3, $4, $4-$3}' > t1
awk '/^t_/ {print $2,$4,$5}' run2/prior_distribution/out | sed "s/,//g;s/)//g;s/ /\t/g" | awk -v OFS="\t" '{print $1, $2, $3, $3-$2}' > t2
paste t* | sed "1i node\tmean_prior_run1\tlower_95%CI_run1\tupper_95%CI_run1\tCI_width_run1\tmean_prior_run2\tlower_95%CI_run2\tupper_95%CI_run2\tCI_width_run2" > prior_times.txt
rm t*

#compare the prior means of each node between two runs
#calculate pearson correlation between two runs
cat prior_times.txt | ${DIR_CSVTK}/csvtk -t corr -T -f mean_prior_run1,mean_prior_run2 1>>Pearson_correlation_coefficient.txt 2>&1
CORR=$(cat Pearson_correlation_coefficient.txt | tail -n 1 | cut -f3)
${DIR_CSVTK}/csvtk -t plot line prior_times.txt -x mean_prior_run1 -y mean_prior_run2 --format pdf -o convergence_plot_prior.pdf

echo -e "Pearson correlation coefficient of the prior means between two runs is "$CORR". You can view the convergence plot 'convergence_plot_prior.pdf' in the folder 3-divergence to check the convergence. The points should fall almost perfectly on the straight line; if not, increase the generation number of MCMC runs."
echo

#generate infinite-sites plot assess whether collecting additional molecular data would improve the analysis
cat prior_times.txt | ${DIR_CSVTK}/csvtk -t corr -T -f mean_prior_run1,CI_width_run1 1>>Pearson_correlation_coefficient.txt 2>&1
${DIR_CSVTK}/csvtk -t plot line prior_times.txt -x mean_prior_run1 -y CI_width_run1 --format pdf -o infinite-sites_plot_run1_prior.pdf
cat prior_times.txt | $DIR_CSVTK/csvtk -t corr -T -f mean_prior_run2,CI_width_run2 1>>Pearson_correlation_coefficient.txt 2>&1
${DIR_CSVTK}/csvtk -t plot line prior_times.txt -x mean_prior_run2 -y CI_width_run2 --format pdf -o infinite-sites_plot_run2_prior.pdf

echo "View the infinite-sites plots of each run for prior distribution. It helps to assess whether collecting additional molecular data would improve the analysis. As the number of sites and loci tend to infinity, a plot of mean prior times vs. credibility interval widths will tend to a straight line. Pearson correlation coefficient of two variables are also calculated and kept in the file 'Pearson_correlation_coefficient.txt'"
echo

#mean posterior vs mean prior
cut -f2,6 posterior_times.txt > t1
cut -f2,6 prior_times.txt > t2
paste t1 t2 > t3
cat t3 | ${DIR_CSVTK}/csvtk -t corr -T -f mean_posterior_run1,mean_prior_run1 1>>Pearson_correlation_coefficient.txt 2>&1
${DIR_CSVTK}/csvtk -t plot line t3 -x mean_posterior_run1 -y mean_prior_run1 --format pdf -o posterior_vs_prior_plot_run1.pdf
cat t3 | ${DIR_CSVTK}/csvtk -t corr -T -f mean_posterior_run2,mean_prior_run2 1>>Pearson_correlation_coefficient.txt 2>&1
${DIR_CSVTK}/csvtk -t plot line t3 -x mean_posterior_run2 -y mean_prior_run2 --format pdf -o posterior_vs_prior_plot_run2.pdf
rm t*

echo "View the mean posterior vs mean prior plots of each run. It helps to assess that whether prior distribution estimated from fossil calibrations are consistent with those posterior dustribution. Pearson correlation coefficient of two variables are also calculated and kept in the file 'Pearson_correlation_coefficient.txt'"

