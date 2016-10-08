#1/bin/bash



for dir in $@
do

    nohup  quast.py -o $dir/both_spades_output_without_error_correction/quast_report $dir/both_spades_output_without_error_correction/contigs.fasta &

done
