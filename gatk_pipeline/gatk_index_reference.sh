DICT=$(echo $1 | sed 's/\..*//' -).dict
bwa index -a bwtsw $1
samtools faidx $1
picard-tools CreateSequenceDictionary R= $1 O= $DICT
