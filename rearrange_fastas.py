#!/usr/bin/python3

import sys
from Bio import SeqIO

label = sys.argv[1]
multi_fasta = sys.argv[2]
if ',' in sys.argv[3]:
    nodes = sys.argv[3].split(',')
else:
    nodes = [sys.argv[3]]

output = open(label + '.fasta', 'w')


def rev_comp(seq):
    """Reverses, complements and returns sequence"""
    rev_seq = seq[::-1]
    compliment_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    rev_comp_seq = ''
    for nuc in rev_seq:
        if nuc in ['A', 'T', 'G', 'C']:
            rev_comp_seq += compliment_dict[nuc]
    return rev_comp_seq

new_sequence = ''
total_length = 0

for node in nodes:
    for seq_record in SeqIO.parse(multi_fasta, "fasta"):
        if node[1:] in str(seq_record.id):
            if node.startswith('-'):
                print (str(seq_record.id) + '--> Reversed')
                sequence = rev_comp(str(seq_record.seq))
            else:
                sequence = str(seq_record.seq)
                print (str(seq_record.id))
            new_sequence+=sequence
            total_length += len(sequence)

print ('new sequence total length = ', total_length)
output.writelines(">" + label + '_' + str(nodes) + '_' + str(total_length) + '_' '\n' + new_sequence)