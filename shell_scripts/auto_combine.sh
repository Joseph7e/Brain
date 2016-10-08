#1/bin/bash

for dir in $@
do
    cat $dir/both_paired-1.fastq.gz >> combined_geminates/atlantic_combined-1.fastq.gz
    cat $dir/both_paired-2.fastq.gz >> combined_geminates/atlantic_combined-2.fastq.gz
	FORWARD=combined_geminates/atlantic_combined-1.fastq.gz;
	REVERSE=combined_geminates/atlantic_combined-2.fastq.gz;
	FORWARD=$(echo $FORWARD);
	REVERSE=$(echo $REVERSE);
        echo $dir
        echo $FORWARD
        echo $REVERSE
done


