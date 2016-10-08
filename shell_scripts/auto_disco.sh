#1/bin/bash

for dir in $@
do
	FORWARD=$dir/both_paired-1.fastq.gz;
	REVERSE=$dir/both_paired-2.fastq.gz;
	FORWARD=$(echo $FORWARD);
	REVERSE=$(echo $REVERSE);

	nohup DiscovarDeNovo NUM_THREADS=1 READS=$FORWARD,$REVERSE OUT_DIR=$dir/disco_assembly &
done
