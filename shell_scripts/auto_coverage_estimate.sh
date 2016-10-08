#1/bin/bash



for dir in Sample_NEM*/
do
    contigs=$dir'both_spades_output_error_correction/contigs.fasta'

    grep ">" $contigs | awk 'BEGIN { FS = "_" } ; {ratio += $4 * $6; total += $4} END { print ratio/total, ratio, total }' >> nemertea_coverage_stats

done
