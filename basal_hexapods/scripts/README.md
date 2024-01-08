## Codes and commands

All commands are labelled as bold. Operating system is popular Linux ones with ‘BASH’ shell, such as Centos, UBUNTU etc.

### 1. Four-cluster likelihood mapping analysis

      iqtree -s XXX.fas -m C60+F+R -lmap ALL -lmclust cluster.nexus -n 0 -T $THREADS

*specify a NEXUS file containing taxon clusters (see **cluster.nexus** file) for quartet mapping analysis*

### 2. Site-wise likelihood analysis

**#calculate site-wise likelihood scores for hypotheses H1, H2 and H3**

      iqtree -s XXX.fas -m EX_EHO+F+R4 -z H1.tre -n 0 -wsl -T $THREADS
      iqtree -s XXX.fas -m EX_EHO+F+R4 -z H2.tre -n 0 -wsl -T $THREADS
      iqtree -s XXX.fas -m EX_EHO+F+R4 -z H3.tre -n 0 -wsl -T $THREADS
      
**#get the number of sites supporting hypotheses H1, or H2, or H3**

      cat H1.sitelh H2.sitelh H3.sitelh | grep "^Tree" | awk '{$1="";print}' | awk 'BEGIN{c=0;} {for(i=1;i<=NF;i++) {num[c,i] = $i;} c++;} END{ for(i=1;i<=NF;i++){str=""; for(j=0;j<NR;j++){ if(j>0){str = str" "} str= str""num[j,i]}printf("%s\n", str)} }' > H1-2-3.sitelh
      awk '{for(i=1;i<=NF;i++)a[NR,i]=$i}END{for(j=1;j<=NF;j++)for(k=1;k<=NR;k++)printf k==NR?a[k,j] RS:a[k,j] FS}' H1-2-3.sitelh > temp1
      bash site-wise_likelihood.sh

### 3. Gene-wise likelihood analysis (i.e., detect distribution of gene tree support)

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
      cat loci.list | grep -v -f loci.incomplete.list > loci.complete.list
      
**#prepare T1.tre, T2.tre and T3.tre from the hypotheses H1, H2 and H3, and prepare all the alignments according to the loci of *loci.complete.list*, and deposited in the *loci/* folder**

      phykit tip_labels T1.tre > species.list
      mkdir compare && cd compare
      bash gene-wise_likelihood.sh

### 4. Phylogenetic inference

All maximum likelihood (ML) supermatrix analyses are performed using IQ-TREE. Mixture model CAT-GTR is performed using PhyloBayes MPI v1.9. Multispecies coalescent (MSC) model is executed using ASTER.

#### 4.1 Partitioned ML model

      iqtree -s XXX.fas -p XXX.partition -m MFP --mset LG --msub nuclear --rclusterf 10 -B 1000 --alrt 1000 -T $THREADS

#### 4.2 PMSF(C60) mixture model

      iqtree -s XXX.fas -m LG+C60+F+R -ft guide.tree -B 1000 -alrt 1000 -T $THREADS

*resulting tree (file **XXX.treefile**) from previous partitioned ML analysis can be used as the guide tree*

#### 4.3 Dayhoff6-recoding model

**#recoding matrix for iqtree and PhyloBayes using Phylogears2**
      
      pgrecodeseq --type=ANY "ARNDCQEGHILKMFPSTWYVX-01223220144145000554?" XXX.fas Dayhoff6.fas
      pgrecodeseq --type=ANY "ARNDCQEGHILKMFPSTWYVX-ABCCDCCABEEBEFAAAFFE?" XXX.fas Dayhoff6.phy

**#IQ-TREE**
      
      iqtree -s Dayhoff6.fas -m GTR+R -B 1000 --alrt 1000 -T $THREADS

#### 4.4 CAT-GTR mixture model

**#transform FASTA format into PHYLIP format using one component of trimAl v1.4.1.**
      
      readal -in XXX.fas -out XXX.phy -phylip_paml

**#PhyloBayes**

      mpirun -np $THREADS pb_mpi -d XXX.phy -cat -gtr -x 1 -1 $CHAIN_NAME
      
**#check convergence and generate the statistics**

      bpcomp -x 8000 1 chain1 chain2
      tracecomp -x 8000 1 chain1 chain2
      
*number '8000' as the burn-in*

#### 4.5 Calculate concordance factors gCF and sCF given a tree 'XXX.treefile'

      iqtree -t XXX.treefile --gcf all.gene.tre --prefix concord
      iqtree -te concord.cf.tree -s XXX.fa --scfl 100 --prefix concord2 -nt $THREADS

#### 4.6 Topology tests

We tested the resultant four alternative topologies with the all four matrices using approximately the unbiased (AU), weighted Kishino-Hasegawa (WKH), and weighted Shimodaira-Hasegawa (WSH) tests under the across-site compositional heterogeneity model in IQ-TREE.

**#merge alternative topologies into a topology hypotheses file 'trees'**
      
      cat *treefile > trees
      
**#topology tests under across-site compositional heterogeneity model**
      
      iqtree -s XXX.fas -m C60+F+R -z trees -n 0 -zb 10000 -zw -au -T $THREADS

### 5. Bayesian leave-one-out cross-validation (LOO-CV) analysis

**#transform FASTA format into PHYLIP format using one component of trimAl v1.4.1.**

      readal -in XXX.fas -out XXX.phy -phylip_paml

**#CAT-GTR model**

      mpirun -np $THREADS pb_mpi -d XXX.phy -cat -gtr -x 1 -1 $CHAIN_NAME

**#LG model**

      mpirun -np $THREADS pb_mpi -d XXX.phy -ncat 1 -lg -x 1 -1 $CHAIN_NAME

**#get the site-specfic log likelihood statistics for all chains**

      mpirun -np $THREADS readpb_mpi -x 1000 10 -sitelogl $CHAIN_NAME     

*numbers ‘1000 10’ as the burn-in of 1000, taking every 10 trees*

**#collect the site-specific scores for all models**

      python3 ~/install/pbmpi-1.9/scripts/read_loocv_waic.py *.sitelogl


## Contact

For questions, please send emails to Prof. Feng Zhang (xtmtd.zf@gmail.com) or Shiyu Du (zjjhdsy@126.com).

