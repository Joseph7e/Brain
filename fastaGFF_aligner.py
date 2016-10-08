#!usr/bin/python3
#Author: Joseph Sevigny
#Affiliation: Hubbard Center For Genome Studies
#Purpose: Check the annotations of MITOS output
#USAGE: python3 .....


import os, sys, subprocess
from Bio import SeqIO

FASTA_dir = sys.argv[1]
GFF_dir = sys.argv[2]
my_genes = sys.argv[3]
flanking = 0
sample_identity = "ANN_"

my_genes_list = my_genes.split(",")
print (my_genes_list)
sys.exit()

output_dir = "alignments/"
single_genes_dir = output_dir + "single_genes"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
if not os.path.exists(single_genes_dir):
    os.mkdir(single_genes_dir)

def extract_file_names(directory,search_term):
    '''extracts a set of file names that contain the
    search term in their name and returns it as a list
    '''
    filenames = []
    file_list = os.listdir(directory)
    for file_name in file_list:
        if search_term in file_name:
            filenames.append(file_name)
    return filenames

gffs = extract_file_names(GFF_dir, "gff")
fastas = extract_file_names(FASTA_dir, "mito")

#create a dictionary of fastas and matching annotations:
fa_gff_dict = {} ####### <-------- The file we will use for alignments and annotations


for fa in fastas:
    for gff in gffs:
        tmp_gff,tmp_ext = gff.split(".")
        if tmp_gff in fa:
            fa_gff_dict[fa] = gff

#Rearrange fasta starting at cox1 and all in the same direction


#concatenate files into a multi_fasta
concat_fasta =  output_dir + sample_identity  + "concat_fastas"
with open(concat_fasta, 'w') as outfile:
    for fname in fastas:
        with open(FASTA_dir + fname) as infile:
            for line in infile:
                outfile.write(line)



def reverse_this(seq):
    """ reverses a sequence and returns it to the user"""
    r_seq = seq[::-1]
    return r_seq

def complement_this(seq):
    """ constructs a complement of the sequence and returns it to the user"""
    compliment_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    rev_seq = ''
    for nuc in seq:
        if nuc in ['A', 'T', 'G', 'C']:
            rev_seq += compliment_dict[nuc]
    return rev_seq

def rev_comp(seq):
    """Combines the reverse_this and complement_this function into one"""
    return(complement_this(reverse_this(seq)))

#Use a gff file instead, it is tab delineated so

def parse_gff(gff_line):
    header, gene_identity, start, stop, direction = '', '', '', '', ''
    try:
        header, program, gene_identity, start, stop, tmp, direction, tmp2, tmp3 = line.split("\t")
    except ValueError:
        header, program, gene_identity, start, stop, tmp, direction, tmp2 = line.split("\t")

    return header, gene_identity, start, stop, direction



#master_file = open("genes_lined_up.joe", 'w')
gene_dict = {}

flanking_dict = {}

print (fa_gff_dict)
for current_fasta, current_gff in fa_gff_dict.items():
    print (current_fasta)
    g = open(GFF_dir + current_gff, 'r')
    for line in g.readlines():
        header, gene_identity, start, stop, direction = parse_gff(line)
        start = int(start)
        stop = int(stop)
        for seq_record in SeqIO.parse(FASTA_dir + current_fasta, "fasta"):
            #add in addition about parsing a list of multiple gene, make a fasta for each gene and a concat of all genes, alignments for all involved
            if my_gene in gene_identity:
                print (my_gene)
                annotation = seq_record.seq[start-1-flanking:stop+flanking]
                annotation = str(annotation)
                left_flank = annotation[:flanking]
                right_flank = annotation[-flanking:]
                flanking_dict[current_fasta] = [left_flank, right_flank]
                if direction == '+':
                    gene_dict[current_fasta] = annotation
                elif direction == '-':
                    annotation = rev_comp(annotation)
                    gene_dict[current_fasta] = annotation


#Create an alignment of fastas_starting at cox1 with

def muscle_alignment(input_multi_fastas, output_name):
    '''
    example_command = ["muscle", "-clw", "-in", "amino_fasta", "-out", "output_a