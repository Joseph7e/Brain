#1/bin/bash

reference=$1
vcf=$2
output=$3

 GATK \
   -T FastaAlternateReferenceMaker \
   -R $reference \
   -o $3 \
   -V $2
