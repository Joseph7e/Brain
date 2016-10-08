# now call haplotypes
ARGS=''
if [ "$#" -gt 2 ]; then
    for i in ${@:2}
    do
        ARGS+='-I '$i' '
    done
else
    ARGS='-I '$2
fi
    #-I ${@:2} \
GATK -T HaplotypeCaller \
    -R $1 \
    $ARGS \
    --genotyping_mode DISCOVERY \
    --sample_ploidy 2 \
    -stand_emit_conf 10 \
    -stand_call_conf 10 \
    --min_base_quality_score 30 \
    --max_alternate_alleles 20 \
    -o raw_variants.vcf 
