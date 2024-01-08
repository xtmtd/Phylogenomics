#!/bin/bash
#2021.11.22 wrote by ZF, 2023.11.14 modified by DSY

#Type 'bash gene-wise_likelihood.sh', e.g. bash gene-wise_likelihood.sh

COMPARE_fun() {
	  mkdir $1 && cd $1
	    cp ../../loci/$1 .
	      cat $1 | grep "^>" | sed "s/^>//g" > taxa.list
	        cat ../../species.list | grep -v -f taxa.list > prune.list
		    phykit prune_tree ../../H1-guide.tree prune.list -o H1.prune.tre
		    phykit prune_tree ../../H2-guide.tree prune.list -o H2.prune.tre
		    phykit prune_tree ../../H3-guide.tree prune.list -o H3.prune.tre
		      cat H1.prune.tre H2.prune.tre H3.prune.tre > candidate.trees
		        iqtree -s $1 -m EX_EHO+F+R4 -z candidate.trees -n 0 -zb 10000 -zw -au -T 1
			  #summary logL and p=value of AU testsi
			    cat $1.iqtree | grep -E ' \+ | - ' | grep "^ " | awk -v OFS='\t' '{print $1,$2,$3,$17}' > $1.summary
			      TREE=$(cat $1.summary | grep -P "\t""0""\t" | cut -f1)
		          echo $1 >> ../loci.T"$TREE".sig
					    cd ..
								}

export -f COMPARE_fun
cat ../loci.complete.list | parallel -I% -j 30 --max-args 1 COMPARE_fun %
