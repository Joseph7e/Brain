#1/bin/bash



for dir in $@
do
    cat $dir/unpaired-* > $dir/unpaired.fastq.gz
    UNPAIRED=$dir/unpaired.fastq.gz
    FORWARD=$dir/paired-1.fastq.gz;
    REVERSE=$dir/paired-2.fastq.gz;
    FORWARD=$(echo $FORWARD);
    REVERSE=$(echo $REVERSE);

    spades.py -1 $FORWARD -2 $REVERSE -s $UNPAIRED -o $dir/spades_output_unpaired

done
