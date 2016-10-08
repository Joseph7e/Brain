#!/usr/bin/python3

import sys, re
from Bio import SeqIO
#Purpose: input two files 1. node_ids\tcoverage\tgc 2. contigs.fasta
#	  output = fasta file with only retained contigs
#    also input gc and coverage cutoffs as sys.argv

table = open(sys.argv[1], 'r')

contig_file = sys.argv[2]

gc_upper = 50#60#32 #50
gc_lower = 37#52#38 #37
coverage_upper = 70#9.5#100#90 #70
coverage_lower = 9#2#10#6 #9

length_cutoff = 0 #500

out_add = 'Nem_b7'

log = open('extracted_nodes_'+out_add+'.txt', 'w')
output_contigs = open('new_contig_file'+out_add+'_'+str(gc_upper)+'_'+str(gc_lower)+'_'+str(coverage_upper)+'_'+str(coverage_lower)+'.fasta', 'w')



good_node_list = []

for line in table.readlines():
    n, length, coverage, gc = line.split('\t')
    node = 'NODE_' + n + '_'
    gc = float(gc.rstrip()); coverage = float(coverage); length = float(length)
    if gc >= gc_lower and gc <= gc_upper and coverage >= coverage_lower and coverage <= coverage_upper and length >= length_cutoff:
        good_node_list.append(node)
        log.writelines(line)

for seq_record in SeqIO.parse(contig_file, "fasta"):
    id = re.findall(r"(NODE_.*)length", seq_record.id)[0]
    if id in good_node_list:
        output_contigs.write(">" + str(seq_record.id) + "\n")
        output_contigs.write(str(seq_record.seq) + "\n")
    else:
        continue


