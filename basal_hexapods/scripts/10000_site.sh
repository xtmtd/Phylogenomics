#!/bin/bash
#2021.11.27 by DSY

#Type 'bash 10000_sites.sh supermatrix.fas', e.g. bash 10000_sites.sh FcC_supermatrix.fas

cp $1 supermatrix.fas

mkdir temp && cd temp
number=$(cat ../supermatrix.fas | wc -l)
for i in `seq 1 2 $number`
do
	j=$(($i+1))
	fn=$(sed -n "$i p" ../supermatrix.fas | sed "s/^>//g")
	seq=$(sed -n "$i p" ../supermatrix.fas | awk -F ">" '{print $2}')
	sed -n "$i,$j p" ../supermatrix.fas > $fn
done
cd ..

cat temp/* > temp1
ls temp/ > species.list
sed -e '/>/d' temp1 > temp2
cat temp2 | sed 's/.\{1\}/& /g' > temp3
cat temp3 | awk '{$1="";print}' | awk 'BEGIN{c=0;} {for(i=1;i<=NF;i++) {num[c,i] = $i;} c++;} END{ for(i=1;i<=NF;i++){str=""; for(j=0;j<NR;j++){ if(j>0){str = str" "} str= str""num[j,i]}printf("%s\n", str)} }' > temp4
shuf -n10000 temp4 > seq_cv

sed 's/[ ][ ]*/,/g' seq_cv > tmp1
awk -F, '{for(i=1;i<=NF;i=i+1){a[NR,i]=$i}}END{for(j=1;j<=NF;j++){str=a[1,j];for(i=2;i<=NR;i++){str=str " " a[i,j]}print str}}' tmp1 > tmp2
sed 's/[ ][ ]*//g' tmp2 > tmp3
rm -rf temp1 temp2 temp3 temp4 temp/ tmp1 tmp2 seq_cv

sed "s/^/>/g" species.list > species1.list
test=$(cat species1.list | wc -l)
for i in `seq 1 1 $test`
do
	j=$(($i+1))
	name=$(sed -n "$i p" species1.list)
	seq=$(sed -n "$i p" tmp3)
	echo $name >> 10000_sites.fas
	echo $seq >> 10000_sites.fas
done
rm tmp3 supermatrix.fas species*
