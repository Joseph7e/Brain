#!/usr/bin/python3

from Bio import SeqIO
import sys


input = sys.argv[1] #"all_16s.fasta"

output = open("new_seqs.fasta", 'w')

cutoff = 1000

for seq_record in SeqIO.parse(input, "fasta"):
    my_seq = (str(seq_record.seq))
    header = (str(seq_record.id))
    if len(my_seq) > cutoff:
        output.writelines("\n>" + header + "\n" + my_seq)