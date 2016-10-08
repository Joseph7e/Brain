#nohup blastn -query $QUERY -subject $dir/spades_output/contigs.fasta -out $dir/spades_output/mito_blast
#!/bin/bash

for dir in $@

do
QUERY=''
    if $dir == SAMPLE_ANN*/
        then
            $QUERY=annelids.fasta
    elif $dir == SAMPLE_NEM*/
        then
            $QUERY=nemertea.fasta
    else
        $QUERY=nematode.fasta
    fi

    echo $QUERY
#nohup blastn -query $QUERY -subject $dir/spades_output/contigs.fasta -out $dir/spades_output/mito_blast 

done
