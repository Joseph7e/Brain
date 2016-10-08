#!/usr/bin/python3
import sys
from Bio import SeqIO

input_handle = open(sys.argv[1], "rU")
output_handle = open(sys.argv[1] + ".fasta", "w")

sequences = SeqIO.parse(input_handle, "genbank")
count = SeqIO.write(sequences, output_handle, "fasta")

output_handle.close()
input_handle.close()
print("Converted %i records" % count)


