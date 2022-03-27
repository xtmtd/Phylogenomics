#this script performs quality control for the raw sequencing data, including removing duplidates and low-quality regions, normalization, error correction.
#Type 'sh quality_control.sh forward_reads_file reverse_reads_file', e.g. sh quality_control.sh illumina.R1.fq.gz illumina.R2.fq.gz


echo "Checking the package dependency......" | tee -a quality_control.log
echo -e "\n" >> quality_control.log


#check pigz
if [ $(which pigz) ]
    then
      echo "pigz ...... OK" | tee -a quality_control.log
      EXE_PIGZ=$(which pigz)
      DIR_PIGZ=${EXE_PIGZ%/*}
    else
      until [ -x $DIR_PIGZ/pigz ]
        do
          read -p "Pigz is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_PIGZ
        done
      echo "pigz ...... OK" | tee -a quality_control.log
fi

#check BBtools
if [ $(which bbduk.sh) ]
    then
      echo "BBtools ...... OK" | tee -a quality_control.log
      EXE_BBTOOLS=$(which bbduk.sh)
      DIR_BBTOOLS=${EXE_BBTOOLS%/*}
    else
      until [ -x $DIR_BBTOOLS/bbduk.sh ]
        do
          read -p "BBtools is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/bbtools):      " DIR_BBTOOLS
        done
      echo "BBtools ...... OK" | tee -a quality_control.log
fi


#group overlapping reads into clumps and remove duplicates
echo "Clumpify and remove duplicates......" | tee -a quality_control.log
$DIR_BBTOOLS/clumpify.sh in1=$1 in2=$2 out1=1.clumped.fq.gz out2=2.clumped.fq.gz pigz dedupe 1>>quality_control.log 2>&1
echo -e "\n" 


#Quality trimming
echo "Quality trimming......" | tee -a quality_control.log
$DIR_BBTOOLS/bbduk.sh in1=1.clumped.fq.gz in2=2.clumped.fq.gz out1=1.trim.fq.gz out2=2.trim.fq.gz ziplevel=5 pigz ordered qtrim=rl trimq=20 minlen=15 ecco=t maxns=5 trimpolya=10 trimpolyg=10 trimpolyc=10 ref="$DIR_BBTOOLS"/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo 1>>quality_control.log 2>&1
echo -e "\n" >> quality_control.log

rm *clumped*

