#1/bin/bash

for dir in $@
do
	cat $dir/*R1*.fastq.gz > $dir/combined_1.fastq.gz
        cat $dir/*R2*.fastq.gz > $dir/combined_2.fastq.gz
        FORWARD=$dir/combined_1.fastq.gz;
	REVERSE=$dir/combined_2.fastq.gz;
	FORWARD=$(echo $FORWARD);
	REVERSE=$(echo $REVERSE);

	echo trimmomatic PE -threads 16 $FORWARD $REVERSE $dir/adapter_paired-1.fastq.gz $dir/_adapter_unpaired-1.fastq.gz $dir/adapter_paired-2.fastq.gz $dir/adapter_unpaired-2.fastq.gz ILLUMINACLIP:/opt/Trimmomatic-0.32/adapters/NexteraPE-PE.fa LEADING:3 TRAILING:3 MINLEN:36
done
