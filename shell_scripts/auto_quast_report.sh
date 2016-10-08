#1/bin/bash


for dir in $@
do
    nohup  quast.py -f -e -o $dir/spades_output_careful/quast_report $dir/spades_output_careful/contigs.fasta &

done
