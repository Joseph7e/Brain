#1/bin/bash

for dir in $@
do
	FORWARD=$dir/all_combined_R1.fastq.gz;
	REVERSE=$dir/all_combined_R2.fastq.gz;
	FORWARD=$(echo $FORWARD);
	REVERSE=$(echo $REVERSE);

	nohup trimmomatic PE -threads 1 $FORWARD $REVERSE $dir/paired-1.fastq.gz $dir/unpaired-1.fastq.gz $dir/paired-2.fastq.gz $dir/unpaired-2.fastq.gz ILLUMINACLIP:/opt/Trimmomatic-0.32/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &
done
