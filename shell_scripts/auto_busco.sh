#1/bin/bash



for dir in $@
do

    #nohup  quast.py -o $dir/both_spades_output_error_correction/quast_report $dir/both_spades_output_error_correction/contigs.fasta &
    nohup python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -l ~/BUSCO_v1.1b1/metazoa -m genome -o busco_output -in $dir/spades_output_careful/contigs.fasta -f &

done
