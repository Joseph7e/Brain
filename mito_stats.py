#!/usr/bin/python3

#purpose: input mito_genome directory, output genome stats

import os, sys
from Bio import SeqIO

input_dir = sys.argv[1]
input_list = os.listdir(input_dir)
output = open(sys.argv[2] + "genome_stats", 'w')



def read_fasta(file_name):
    """reads a fasta file and returns a dictionary with headers as keys and sequences as values"""
    cfo = open(file_name, 'r')
    fasta_dict = {}
    for line in cfo.readlines():
        if line[0] == ">":
            header = line
        else:
            if header not in fasta_dict:
                fasta_dict[header] = line.rstrip()
            else:
                fasta_dict[header] += line.rstrip()
    return fasta_dict


print("Creating a dictionary of contig headers and sequences")

output.writelines("#USE --> sort -t $'TAB' -nr -k3 details.file > sorted.file , To sort the file by the specified column (replace '3' with '2' or '4' for others)\n#\n#NODE_NAME\tfile\tLength\tGC_Percentage\n")

for file in input_list:
    if ".fa" in file or "contig" in file:
        file = input_dir+file
        avg_list = []
        all_gc = 0
        length = 0
        contig_dict = read_fasta(file) #Uses the read_fasta definition to create a LARGE contig dictionary file
        for head, seq in sorted(contig_dict.items()): #sorts the contig nodes alphabetically
            node_length = len(seq)
            node = head

            length += node_length

            G = seq.count('G'); C = seq.count('C') #creates a count for each nucleotide in each of the expanded sequences
            GC_content = ((G+C)/len(seq)*100)
            avg_list.append(GC_content)
            output.writelines("{0}\t{1}\t{2:.2f}\n".format(file,node_length,GC_content))
            all_gc += GC_content
        avg_gc = all_gc/len(avg_list)
        output.writelines(file+"\t"+str(length)+"\t"+str(avg_gc)+"\n")

