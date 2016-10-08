# get sample name for adding read groups
SAMPLE_NAME=$(echo $2 | sift --replace '$1' '(([^_]*_){3})([^_]+)' | sed 's/.$//' -)
BARCODE=$(echo $2 | sift --replace '$1' '.*_([^-]*-[^_]*)' -)

# run bwa
bwa mem -M -t 8  $1 $2 $3 > $4
BASE_OUT=$(echo $4 | sed "s/\..*//" -)

# pipe into sam processing
samtools view -@ 8 -Sb -F 4  $4  | samtools sort -@ 8 - $BASE_OUT
samtools index $BASE_OUT.bam


# begin marking and adding read groups
#picard-tools MarkDuplicates I= $BASE_OUT.bam O= marked.bam M= "$BASE_OUT"_metrics.txt
mv $BASE_OUT.bam marked.bam
picard-tools AddOrReplaceReadGroups \
    RGLB= illumLib \
    RGPU= $BARCODE \
    RGSM= $SAMPLE_NAME \
    RGPL= illumina   \
    I= marked.bam \
    O= grouped.bam

# remove intermediates
rm $4
rm $BASE_OUT.bam
rm marked.bam
mv grouped.bam $BASE_OUT.bam
samtools index $BASE_OUT.bam

# now call haplotypes
