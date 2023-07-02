#!/bin/bash
#2021.12.20 by DSY

#Type 'bash site-wise_likelihood.sh', e.g. bash site-wise_likelihood.sh

#the number '201896' represents the total number of sites for this matrix, users should modify this number according to the own dataset (command to check this number: cat H1-2-3-4.sitelh | wc -l)
mkdir likelihood
for i in `seq 1 1 201896`
 do
   awk '{print $1}' temp$i > likelihood/$i.fas
   cat likelihood/$i.fas | cat -n | sort -k2 | head -1 > likelihood/$i.sort.fas
   rm likelihood/$i.fas
   j=$(($i+1))
   awk '{$1="";print $0}' temp$i > temp$j
   rm temp$i
 done

ls likelihood/ > list
for list in $(cat list)
 do 
  cat likelihood/$list >> likelihood.list
 done

awk '{print $1}' likelihood.list > number
cat number | grep "1" | wc -l > 1
cat number | grep "2" | wc -l > 2
cat number | grep "3" | wc -l > 3
cat number | grep "4" | wc -l > 4
cat 1 2 3 4 | cat -n > site-wise

rm -rf list temp* likelihood/ number 1 2 3 4
