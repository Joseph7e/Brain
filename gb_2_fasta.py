#!/usr/bin/python3
from Bio import SeqIO
import re, sys

converted_output_straight_fromBio = open('trainingset.fasta', 'w')

def gb2fasta(gb, fasta_name):
    input = open(gb, "rU")
    output = open(fasta_name, "w")
    sequences = SeqIO.parse(input, "genbank")
    count = SeqIO.write(sequences, output, "fasta")
    input.close()
    output.close()
    print ("Converted {} record(s)".format(count))
    converted_output_straight_fromBio.writelines(output)
    return output

gb2fasta(sys.argv[1], "output.fasta")

def get_gene_location(gb):
    gene_start = 0
    gene_stop = 0
    for line in open(gb, 'r'):
        if "CDS" in line:
            line = line.replace("..", ":")
            line = line.replace(")", ":")
            line = line.replace("(", ":")
            tmp1, start, stop, tmp2 = line.split(":")
            gene_start = int(start)
            gene_stop = int(stop)
    return gene_start, gene_stop



gene_start, gene_stop = get_gene_location(sys.argv[1])

print (gene_start, gene_stop)

def extract_region(fasta, start, stop, output_name):
    with open(output_name, "w") as my_fasta:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            my_fasta.write(">" + str(seq_record.id) + "\n")
            my_fasta.write(str(seq_record.seq[start-1:stop-1] + "\n"))



extract_region("output.fasta", gene_start, gene_stop, "new_fasta")

