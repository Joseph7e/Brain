#1/bin/bash



for dir in $@
do
    FORWARD=$dir/*aired-1*.fastq.gz;
    REVERSE=$dir/*aired-2*.fastq.gz;
    FORWARD=$(echo $FORWARD);
    REVERSE=$(echo $REVERSE);

    nohup ~/SPAdes-3.6.0-Linux/bin/spades.py --only-assembler -1 $FORWARD -2 $REVERSE -o spades_assembly &

done
