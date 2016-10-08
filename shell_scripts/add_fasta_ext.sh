#1/bin/bash



for fasta in ~/MITOCHONDRIAL_GENOMES/annelids/*contig*
do
    echo $fasta
    mv $fasta $fasta.fasta
done

