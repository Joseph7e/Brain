#!/usr/bin/python3

import sys, re
from Bio import SeqIO

input_fasta = sys.argv[1]

for seq_record in SeqIO.parse(open(input_fasta), "fasta"):
    header = ">" + input_fasta.replace('.fas','')[-2:] + "_" +seq_record.id
    my_seq = str(seq_record.seq)
    print (seq_record.id)
    gene = re.findall(r"[\+-];,(.*)", seq_record.id)
    print (gene)
    output = open(gene[0] + ".fas", 'a')
    output.writelines((header + '\n' + my_seq + '\n'))