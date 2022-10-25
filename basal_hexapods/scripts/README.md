
## Codes and commands

All commands are labelled as bold. Operating system is popular Linux ones with ‘BASH’ shell, such as Centos, UBUNTU etc.

### 1. Identify 'inconsistent' genes resulting topological differences in supermatrix-based (T1) and supertree-based (T2) phylogenies

**#estimate site-wise loglk for alternative hypothesis**

      mkdir GLS && cd GLS
      cat T1.tre T2.tre > ML2ASTRAL.tre
      iqtree -s FcC_supermatrix.fas -Q FcC_supermatrix_partition.txt -m LG+F+R4 -z ML2ASTRAL.tre -wsl -T 12 --prefix ML2ASTRAL
      perl GLS_parser_v1.pl ML2ASTRAL.sitelh FcC_supermatrix_partition.txt GLS

  *the perl script of **GLS_parser_v1.pl** from Shen et al. (2021)*

**#estimate the quartet score for alternative hypothesis, and prune tips of T1 and T2 to make their tips equal to gene trees**

      mkdir GQS && cd GQS
      GQS_fun() {
        phykit tl trees/$1 > tip.$1
        phykit tl T1.tre > tip.T1.$1
        cat tip.T1.$1 | grep -v -f tip.$1 > tip.missing.$1
        phykit prune T1.tre tip.missing.$1 -o T1.prune.tre.$1
        phykit prune T2.tre tip.missing.$1 -o T2.prune.tre.$1
        java -jar $ASTRAL_PATH/astral.5.7.1.jar -i trees/$1 -q T1.prune.tre.$1 2> log.txt.T1.$1
        GQS_T1=$(cat log.txt.T1.$1 | grep "^Final quartet" | cut -d ":" -f2 | sed "s/ //g")
        java -jar $ASTRAL_PATH/astral.5.7.1.jar -i trees/$1 -q T2.prune.tre.$1 2> log.txt.T2.$1
        GQS_T2=$(cat log.txt.T2.$1 | grep "^Final quartet" | cut -d ":" -f2 | sed "s/ //g")
        diff=$(echo "scale=0;($GQS_T1-$GQS_T2)"| bc)
        echo -e $1"\t"$GQS_T1"\t"$GQS_T2"\t"$diff > trees/$1.GQS
        rm tip*$1 log.txt*$1 *prune.tre.$1
      }
      export -f GQS_fun
      cat tree.list | parallel -I% -j $THREADS --max-args 1 GQS_fun %
      cat trees/*GQS > GQS_table.txt

**#get the inconsistent loci**

      cat GQS/GQS_table.txt | cut -f1 | cut -d "." -f1 > loci.list
      for loci in $(cat loci.list)
        do  
          num1=$(cat GLS/GLS_table.txt | grep $loci | cut -f5 | sed "s/e-//g;s/e+//g")
          num2=$(cat GQS/GQS_table.txt | grep $loci | cut -f4)
          c_num1=`echo "$num1 > 0" |bc`
          c_num2=`echo "$num2 > 0" |bc`
          test "$c_num1" = 1 -a "$c_num2" = 1 && echo $loci >> loci.consistent
          c_num3=`echo "$num1 < 0" |bc`
          c_num4=`echo "$num2 < 0" |bc`
          test "$c_num3" = 1 -a "$c_num4" = 1 && echo $loci >> loci.consistent
        done

### 2. Four-cluster likelihood mapping analysis

      iqtree -s FcC_supermatrix.fas -m LG+C60+F+R -ft guide.tree -lmap ALL -lmclust cluster.nexus -n 0 -T $THREADS

*specify a NEXUS file containing taxon clusters (see **cluster.nexus** file) for quartet mapping analysis*

### 3. Site-wise likelihood analysis

**#calculate site-wise likelihood scores for hypotheses H1, H2, H3 and H4**

      iqtree -s FcC_supermatrix.fas -m EX_EHO+F+R4 -z H1.tre -n 0 -wsl -T $THREADS
      iqtree -s FcC_supermatrix.fas -m EX_EHO+F+R4 -z H2.tre -n 0 -wsl -T $THREADS
      iqtree -s FcC_supermatrix.fas -m EX_EHO+F+R4 -z H3.tre -n 0 -wsl -T $THREADS
      iqtree -s FcC_supermatrix.fas -m EX_EHO+F+R4 -z H4.tre -n 0 -wsl -T $THREADS
      
**#get the number of sites supporting hypotheses H1, or H2, or H3, or H4**

      cat H1.sitelh H2.sitelh H3.sitelh H4.sitelh | grep "^Tree" | awk '{$1="";print}' | awk 'BEGIN{c=0;} {for(i=1;i<=NF;i++) {num[c,i] = $i;} c++;} END{ for(i=1;i<=NF;i++){str=""; for(j=0;j<NR;j++){ if(j>0){str = str" "} str= str""num[j,i]}printf("%s\n", str)} }' > H1-2-3-4.sitelh
      awk '{for(i=1;i<=NF;i++)a[NR,i]=$i}END{for(j=1;j<=NF;j++)for(k=1;k<=NR;k++)printf k==NR?a[k,j] RS:a[k,j] FS}' H1-2-3-4.sitelh > temp1
      bash site-wise_likelihood.sh

### 4. Gene-wise likelihood analysis (i.e., detect distribution of gene tree support)

**#exclude loci with one of the four group taxa (P: Protura; D: Diplura; C: Collembola; I: Insecta) complete lacking**

      for loci in $(cat loci.list)
        do
          seqkit seq -n loci/$loci > taxa
          grep -f P taxa > t1
          grep -f D taxa > t2
          grep -f C taxa > t3
          grep -f I taxa > t4
          test -s t1 && test -s t2 && test -s t3 && test -s t4 && cp loci/$loci $loci || echo $loci >> loci.incomplete.list
          rm t*
        done

**#prepare T1.tre, T2.tre, T3.tre and T4.tre from the hypotheses H1, H2, H3 and H4**

      phykit tip_labels T1.tre > species.list
      mkdir compare && cd compare
      bash gene-wise_likelihood.sh

**#too few significant genes, thus restrict genes of |logL1-logL2| >=2 as strong evidence**

      for loci in $(cat loci.list)
        do
          TREE=$(cat $loci/*.summary | grep -P "\t"0"\t" | cut -f1)
          diff=$(cat $loci/*.summary | grep -v -P "\t"0"\t" | cut -f3 | awk '{printf("%f",$0)}')
          num=`echo "$diff < 2" |bc`
          test "$num" = 0 && echo $loci >> loci.T"$TREE".strong || echo $loci >> loci.T"$TREE".weak
        done

### 5. Phylogenetic inference

All maximum likelihood (ML) supermatrix analyses are performed using IQ-TREE. Mixture model CAT-GTR is performed using PhyloBayes MPI v1.8b. Multispecies coalescent (MSC) model is executed using ASTRAL-III v5.6.1.

#### 5.1 Partitioned ML model

      iqtree -s FcC_supermatrix.fas -p FcC_supermatrix_partition.txt -m MFP --mset LG --msub nuclear --rclusterf 10 -B 1000 --alrt 1000 -T $THREADS

#### 5.2 Across-site compositional heterogeneity model
      
      iqtree -s FcC_supermatrix.fas -m C60+F+R -B 1000 -alrt 1000 -T $THREADS

#### 5.3 PMSF(C60) mixture model

      iqtree -s FcC_supermatrix.fas -m LG+C60+F+R -ft guide.tree -B 1000 -alrt 1000 -T $THREADS

*resulting tree (file **XXX.treefile**) from previous partitioned ML analysis can be used as the guide tree*

#### 5.4 CAT-GTR mixture model

**#transform FASTA format into PHYLIP format using one component of trimAl v1.4.1.**
      
      readal -in FcC_supermatrix.fas -out FcC_supermatrix.phy -phylip_paml

**#PhyloBayes**

      mpirun -np $THREADS pb_mpi -d FcC_supermatrix.phy -cat -gtr -x 1 -1 $CHAIN_NAME
      
**#check convergence and generate the statistics**

      bpcomp -x 3000 1 chain1 chain2
      tracecomp -x 3000 1 chain1 chain2
      
*number '3000' as the burn-in*

#### 5.5 Calculate concordance factors gCF and sCF given a tree 'XXX.treefile'

      iqtree -t XXX.treefile --gcf all.gene.tre -s FcC_supermatrix.fas --scf 100 --prefix concord -T $THREADS

#### 5.6 Topology tests

We tested the resultant four alternative topologies with the all four matrices using approximately the unbiased (AU), weighted Kishino-Hasegawa (WKH), and weighted Shimodaira-Hasegawa (WSH) tests under the across-site compositional heterogeneity model and PMSF(C60) model in IQ-TREE.

**#merge alternative topologies into a topology hypotheses file 'trees'**
      
      cat *treefile > trees
      
**#topology tests under across-site compositional heterogeneity model**
      
      iqtree -s FcC_supermatrix.fas -m C60+F+R -z trees -n 0 -zb 10000 -zw -au -T $THREADS
      
**#topology tests under PMSF(C60) model**
      
      iqtree -s FcC_supermatrix.fas -m LG+C60+F+R -ft guide.tree -z trees -n 0 -zb 10000 -zw -au -T $THREADS

### 6. Bayesian cross-validation analysis

**#select subsample of 10,000 sites for cross-validation analysis**



**#transform FASTA format into PHYLIP format using one component of trimAl v1.4.1.**

      readal -in learn.fas -out learn.phy -phylip_paml
      readal -in test.fas -out test.phy -phylip_paml

**#CAT-GTR model**

      mpirun -np $THREADS pb_mpi -d learn.phy -cat -gtr -x 1 -1 $CHAIN_NAME

**#LG model**

      mpirun -np $THREADS pb_mpi -d learn.phy -ncat 1 -lg -x 1 -1 $CHAIN_NAME

**#**

      mpirun -np $THREADS readpb_mpi -cv test.phy -x 2000 1 5000 $CHAIN_NAME

*numbers ‘2000 1 5000’ as the burn-in of 2000, taking every tree, up to the 5000th point of the chains*

### 7. Phylogeny without outgroup taxa

**#an unrooted tree was inferred by using reversible models**

      iqtree -s without_outgroup.fas -p without_outgroup_partition.txt -B 1000 -T $THREADS

**#a rooted tree with linked non-reversible models was inferred**

      iqtree -s without_outgroup.fas -p unroot.best_model.nex --model-joint NONREV -B 1000 -T $THREADS

**#a rooted tree with linked non-reversible models was inferred to measure the confidence in the root placement**

      iqtree -s without_outgroup.fas -p unroot.best_model.nex --model-joint NONREV --root-test -zb 10000 -zw -au -te rooted_tree.treefile -T $THREADS


## Contact

For questions, please contact fzhang@njau.edu.cn (Prof. Feng Zhang) or zjjhdsy@126.com (Shiyu Du).
