#!/bin/bash

if [ $# -lt 2 ]
then
	echo "Usage: ./generate_alignment.sh <ref_sequence> <out_root> "
	echo "[note]: out_root shall contain aligned signals in '.reso_new' file "
	exit
fi


#========= default parameters ===========#
#========= default parameters ===========#
#-> main
SAMPLE_RATE=8
SAMPLE_LEN=50

#------ read input arguments ---------#
ref_sequence=$1
out_root=$2

#------ mkdir TMP directory ------#
relnam=`basename $1`
txp="${relnam}_PLOT_SIGNAL"
mkdir -p $txp


#------ step 0: pre-process ------#
echo "step 0: pre-process"
python util/genome_preprocess.py -i $ref_sequence -o $txp/refseq.fasta -r 1
ref_len=`tail -n1 $txp/refseq.fasta | wc | awk '{print $3-1}'`
sample_len=`head -n1 $txp/refseq.fasta | awk '{a=int(a*b);print a}' a=$ref_len b=$SAMPLE_LEN`


#------ step 2: transfer ref and snp sequence to signal ------#
echo "step 2: transfer ref and snp sequence to signal"

#-> 2.1 ref sequence
util/seq2sig -i $txp/refseq.fasta -o $txp/ref.seqsig -k 1 -z 1 -N 1
len=`wc $txp/ref.seqsig | awk '{print $1}'`
rm -f $txp/ref.position
for ((i=1;i<=$len;i++))
do
	echo "$i" >> $txp/ref.position
done
paste -d ' ' $txp/ref.seqsig $txp/ref.position > $txp/ref.seqsig_posi
util/NanoRaw_Label $txp/ref.seqsig_posi 1> $txp/ws1 2> $txp/ws2
util/Signal_Transform $txp/ws2 $txp/ws1 $SAMPLE_LEN $txp/refout
rm -f $txp/ws1 $txp/ws2

#-> 2.3 filter candidate_list to generate final_list
rm -f $txp/reso.signal
ls $out_root/*.reso_new | awk -F"/" '{print $NF}' | cut -d '.' -f 1 > $txp/proc_list
for i in `cat $txp/proc_list`
do
	grep -v "#" $out_root/$i.reso_new | awk '{if(NF==6){$1=$1+1; $2=$2+1; print $0}}' > $txp/$i.reso
	awk '{print $5" "$0" "$1}' $txp/$i.reso > $txp/$i.reso_posi
	util/NanoRaw_Label $txp/$i.reso_posi 1> $txp/ws1 2> $txp/ws2
	util/Signal_Transform $txp/ws2 $txp/ws1 $SAMPLE_LEN $txp/tttout
	cat $txp/tttout.signal >> $txp/reso.signal
	rm -f $txp/ws1 $txp/ws2
	rm -f $txp/tttout.*
done

#------ step 4: SNP visualization ---------#
echo "step 4: SNP visualization"
cat $txp/refout.signal $txp/refout.signal $txp/reso.signal > $txp/final_signal
python util/PlotSignal_SNP.py $txp/final_signal $txp/refout.label $relnam.pdf 0 $sample_len

#------ remove temporary files ----#
rm -rf $txp

#=========== exit ==============#
exit 0



