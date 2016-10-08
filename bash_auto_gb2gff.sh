#1/bin/bash



for gb in genbank/*
do
    genbank_to_gff.py $gb
done
