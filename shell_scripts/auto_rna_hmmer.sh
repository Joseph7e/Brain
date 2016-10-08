#1/bin/bash

for fasta in $@
do
    FASTA_OUTPUT=$(basename $fasta).rnammer.fasta;
    GFF_OUTPUT=$(basename $fasta).gff;
    XML_OUTPUT=$(basename $fasta).xml;
    echo $FASTA_OUTPUT
    echo $GFF_OUTPUT
    echo $XML_OUTPUT
    #perl /genome/joseph7e/program_rnammer/rnammer -S bacterial -m tsu,ssu,lsu -gff $GFF_OUTPUT -f $FASTA_OUTPUT -xml $XML_OUTPUT $fasta
    perl /genome/joseph7e/program_rnammer/rnammer -S eukaryotic -m tsu,ssu,lsu -gff $GFF_OUTPUT -f $FASTA_OUTPUT -xml $XML_OUTPUT $fasta
done

