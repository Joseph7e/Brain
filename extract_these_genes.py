#!/usr/bin/python3

#purpose = extract a set of genes from a directory based on line search term

import gzip, os, sys

FASTA_dir = sys.argv[1]
GFF_dir = sys.argv[2]
output_dir = sys.argv[3]

def extract_file_names(search_term, path_to_dir):
    '''extracts a set of file names that contain a search term and returns it as a list'''
    filenames = []
    file_list = os.listdir(path_to_dir)
    for file_name in file_list:
        if search_term in file_name:
           filenames.append(file_name)
    return filenames


def rev_comp(seq):
    """Reverses, complements and returns sequence"""
    rev_seq = seq[::-1]
    compliment_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    rev_comp_seq = ''
    for nuc in rev_seq:
        if nuc in ['A', 'T', 'G', 'C']:
            rev_comp_seq += compliment_dict[nuc]
    return rev_comp_seq



def extract_region(fasta, header_search, start, stop, output_name, direction, header_addition):
    """ PURPOSE: Extracts a specified region from a fasta using header_search, coordinates, and direction(+-),
                        and writes to append able output
        USAGE: extract_region(fasta, header_search, start, stop, output_name, direction, header_addition)"""
    with open(output_name, "a") as my_fasta:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            if header_search in str(seq_record.id):
                my_fasta.write(">" + header_addition + str(seq_record.id) + "\n")
                my_seq = (str(seq_record.seq[start-1:stop]))
                if direction == "-":
                    rev_seq = rev_comp(my_seq)
                    my_fasta.write(rev_seq + "\n")
                else:
                    my_fasta.write(my_seq + "\n")
            else:
                my_fasta.write(">" + header_addition + str(seq_record.id) + "\n")
                my_seq = (str(seq_record.seq[start-1:stop]))
                if direction == "-":
                    rev_seq = rev_comp(my_seq)
                    my_fasta.write(rev_seq + "\n")
                else:
                    my_fasta.write(my_seq + "\n")

def parse_gff(gff_line):
    '''input a line from a gff and return the important attributes'''
    try:
        header, program, gene_identity, start, stop, tmp, direction, tmp2 = gff_line.split("\t")
    except ValueError:
        header, program, gene_identity, start, stop, tmp, direction, tmp2, gene_id = gff_line.split("\t")
        #gene_identity = gene_id[5:11]
    ###try and convert all to the same format
    #gene_identity = gene_identity.lower(); gene_identity = gene_identity.replace("nd", "nad").replace("cyt","co").replace("\n", "")
    return header, gene_identity, start, stop, direction


gffs = extract_file_names("gff", GFF_dir)
fastas = extract_file_names(".fa", FASTA_dir)


fa_gff_dict = {} ####### <-------- The file we will use for alignments and annotations


for fa in fastas:
    for gff in gffs:
        tmp_gff,tmp_ext = gff.split(".")
        tmp_gff = tmp_gff.replace("_", "")
        tmp_gff = tmp_gff.replace("-", "")
        fatmp = fa.replace("_", "")
        fatmp = fatmp.replace("-", "")
        if tmp_gff in fatmp:
            fa_gff_dict[fa] = gff #fills the dictionary with mito fasta and matching gff

print (fa_gff_dict)

sys.exit()
for current_fasta, current_gff in fa_gff_dict.items():
    g = open(GFF_dir + current_gff, 'r')
    print (mess_with_font.Green + "\nDetails for each gene in the file\n" + mess_with_font.ENDC )
    for line in g.readlines():
        #if "trn" in line: #skips if line is not a protein coding gene
        #    continue
        if line[0] == "#":
            continue
        elif "CDS" in line:
            continue
        else:
            header, gene_identity, start, stop, direction = parse_gff(line) #grabs relavent details from gff
            print(header, gene_identity, start, stop, direction)
            start = int(start); stop = int(stop)
        for current_gene in gene_names: #extracts pcg and writes to appendable file
            if current_gene in gene_identity:
                gene_file = output_dir + current_gene + ".fasta"
                extract_region(FASTA_dir + current_fasta, "e", start-flanking, stop+flanking, gene_file, direction, gene_identity + current_fasta)

