#!/usr/bin/python3
import sys
from Bio import SeqIO


input_fasta = sys.argv[1]
output_dir = sys.argv[2]

count = 0

for seq_record in SeqIO.parse(input_fasta, "fasta"):
    count += 1
    add = output_dir[:-1]
    output= open(output_dir + 'gene_' + add + "_" + str(count) + '.faa', 'w')
    header = seq_record.id
    my_seq = str(seq_record.seq)
    output.write(">" + header +"\n" + my_seq)

