#!/usr/bin/python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq
import Bio.Alphabet

fasta = sys.argv[1]
output_name = sys.argv[2]
genetic_code = sys.argv[3] #5

for seq_record in SeqIO.parse(fasta, "fasta"):
    header = ">" + str(seq_record.id) + "\n"
    coding_dna = Seq(str(seq_record.seq))
    protein_seq = coding_dna.translate(table=genetic_code)
    with open(output_name, 'a') as o:
        o.writelines(header)
        o.writelines(protein_seq + "\n")