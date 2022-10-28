#!/bin/bash
#2021.11.22 wrote by ZF, 2021.12.20 modified by DSY

#Type 'bash gene-wise_likelihood.sh', e.g. bash gene-wise_likelihood.sh

COMPARE_fun() {
	  mkdir $1 && cd $1
	    cp ../../loci/$1 .
	      cat $1 | grep "^>" | sed "s/^>//g" > taxa.list
	        cat ../../species.list | grep -v -f taxa.list > prune.list
		    phykit prune_tree ../../T1.tre prune.list -o T1.prune.tre
		    phykit prune_tree ../../T2.tre prune.list -o T2.prune.tre
		    phykit prune_tree ../../T3.tre prune.list -o T3.prune.tre
		    phykit prune_tree ../../T4.tre prune.list -o T4.prune.tre
		      cat T1.prune.tre T2.prune.tre T3.prune.tre T4.prune.tre > candidate.trees
		        iqtree -s $1 -m EX_EHO+F+R4 -z candidate.trees -n 0 -zb 10000 -zw -au -T 1
			  #summary logL and p=value of AU testsi
			    cat $1.iqtree | grep -E ' \+ | - ' | grep "^ " | awk -v OFS='\t' '{print $1,$2,$3,$17}' > $1.summary
			      sig=$(cat $1.summary | grep -P "\t""-$")
			        if [ "$sig" ] ; then
					    TREE=$(cat $1.summary | grep -P "\t""\+$" | cut -f1)
					        echo $1 >> ../loci.T"$TREE".sig
						  else
							      TREE=$(cat $1.summary | grep -P "\t"0"\t" | cut -f1)
							          echo $1 >> ../loci.T"$TREE".exclude_sig
								    fi
								      rm $1 $1.ckp.gz
								        cd ..
								}

export -f COMPARE_fun
cat ../loci.complete.list | parallel -I% -j 30 --max-args 1 COMPARE_fun %
