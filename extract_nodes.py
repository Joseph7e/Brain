#!/usr/bin/python3

# takes a list of node_names and extracts them to multifasta

import sys
from Bio import SeqIO

specimen = sys.argv[1] #output_name
contigs = sys.argv[2]  #contigFile
nodes = sys.argv[3] #   NODE_1_,NODE_359393_

###CONFIGS###

node_list = nodes.split(",")
output_file = specimen + "nodes.fasta"

count = len(node_list)

with open (output_file, 'w') as cfo:
    for current_node in node_list:
        print (current_node)
        if count:
            for seq_record in SeqIO.parse(contigs, "fasta"):
                if current_node in str(seq_record.id):
                    print ("BOOM")
                    count -= 1
                    cfo.write(">" + specimen + current_node + str(seq_record.id) + "\n")
                    cfo.write(str(seq_record.seq) + "\n")
