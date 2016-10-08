#!/usr/bin/python3
import sys
input_gff = "../files/gffs/ANN-B4_mito_contigs.gff"
genbank_fasta = "ANN-B4.genbank"
output = "ANN-B4_mitochondria.gb"

organism = "Annelida Phyllodocidae"
length = "15521"
sample = "Nemertea-A02"




my_string = ''
my_string += 'LOCUS       {}      {} bp    DNA     circular   UNC 13-JAN-2016\nDefinition  {} mitochondrion, complete genome.\nSOURCE      mitochondrion {}\n  ORGANISM  {}\n            Eukaryota; Metazoa; Lophotrochozoa; Nemertea; Enopla;\n            Hoplonemertea; Monostilifera; Ototyphlonemertidae;\n            Ototyphlonemertes\nREFERENCE   1  (bases 1 to {})\nFEATURES             Location/Qualifiers'.format(sample,length,organism,organism,organism,length)

my_string += '\n     source          1..{}\n                     /organism="{}"\n                     /organelle="mitochondrion"\n                     /mol_type="genomic DNA"\n'.format(length,organism)


def parse_gff(gff_line):
    '''input a line from a gff and return the important attributes'''
    try:
        header, program, gene_identity, start, stop, tmp, direction, tmp2 = gff_line.split("\t")
    except ValueError:
        header, program, gene_identity, start, stop, tmp, direction, tmp2, tmp3 = gff_line.split("\t")
    return header, gene_identity, start, stop, direction

    # try:
    #     header, program, gene_id, start, stop, tmp, direction, tmp2, gene_identity, tmp = gff_line.split("\t")
    #     gene_identity = gene_identitytmp[5:11]
    # except ValueError:
    #     header, program, gene_identity, start, stop, tmp, direction, tmp2 = gff_line.split("\t")
    # gene_identity = gene_identity.lower(); gene_identity = gene_identity.replace("nd", "nad").replace("cyt","co").replace("\n", "")
    # return header, gene_identity, start, stop, direction

def construct_annotation(gene_id, direction, start, stop):
    print (gene_id, direction)
    if "trn" in gene_id:
        type = "tRNA"
        if direction == "+":
            line = '     {}            {}..{}\n                     /locus_tag="{}"\n                     /product="{}"\n'.format(type,start, stop, gene_id, gene_id)
        if direction == "-":
            line = '     {}            complement({}..{})\n                     /locus_tag="{}"\n                     /product="{}"\n'.format(type,start, stop, gene_id, gene_id)
    elif "rrn" in gene_id:
        type = "rRNA"
        if direction == "+":
            line = '     {}            {}..{}\n                     /locus_tag="{}"\n                     /product="{}"\n'.format(type,start, stop, gene_id, gene_id)
        if direction == "-":
            line = '     {}            complement({}..{})\n                     /locus_tag="{}"\n                     /product="{}"\n'.format(type,start, stop, gene_id, gene_id)
    else:
        type = "gene"
        if direction == "+":
            line = '     {}            {}..{}\n                     /gene="{}"\n                     /locus_tag="{}"\n'.format(type,start, stop, gene_id, gene_id)
        if direction == "-":
            line = '     {}            complement({}..{})\n                     /gene="{}"\n                     /locus_tag="{}"\n'.format(type,start, stop, gene_id, gene_id)
    return line
for line in open(input_gff, 'r'):
    header, locus_tag, start, stop, direction = parse_gff(line)
    new_line = construct_annotation(locus_tag, direction, start, stop)
    my_string += new_line

print (my_string)

with open(genbank_fasta, 'r') as f:
    flag = False
    for local in f:
        if "ORIGIN" in local:
            flag = True
        if flag:
            my_string += local + "\n"

with open(output, 'w') as file:
    file.writelines(my_string)
