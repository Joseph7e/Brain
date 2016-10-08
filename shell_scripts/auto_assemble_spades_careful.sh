#1/bin/bash

for dir in $@
do
    cat $dir/both_unpaired-1.fastq.gz $dir/unpaired-2.fastq.gz > $dir/unpaired.fastq.gz
    FORWARD=$dir/both_paired-1.fastq.gz;
    REVERSE=$dir/both_paired-2.fastq.gz;
    UNPAIRED=$dir/unpaired.fastq.gz;
    FORWARD=$(echo $FORWARD);
    REVERSE=$(echo $REVERSE);

    spades.py --careful -1 $FORWARD -2 $REVERSE -s $UNPAIRED -o $dir/spades_output_careful_unpaired
    quast.py -o $dir/quast_output_contigs -f --eukaryotes $dir/spades_output_careful_unpaired/contigs.fasta
    BUSCO_v1.2.py -in $dir/spades_output_careful_unpaired/contigs.fasta -l /opt/BUSCO_v1.2/metazoa/ -o $dir/busco_metazoa_long --long
    
done
