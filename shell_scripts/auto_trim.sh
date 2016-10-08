#1/bin/bash

for dir in $@
do
	cat $dir/*R1*.fastq.gz > $dir/combined_1.fastq.gz
        cat $dir/*R2*.fastq.gz > $dir/combined_2.fastq.gz
        FORWARD=$dir/combined_1.fastq.gz;
	REVERSE=$dir/combined_2.fastq.gz;
	FORWARD=$(echo $FORWARD);
	REVERSE=$(echo $REVERSE);

	nohup trimmomatic PE -threads 16 $FORWARD $REVERSE $dir/new_paired-1.fastq.gz $dir/unpaired-1.fastq.gz $dir/new_paired-2.fastq.gz $dir/new_unpaired-2.fastq.gz ILLUMINACLIP:/opt/Trimmomatic-0.32/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &
done
