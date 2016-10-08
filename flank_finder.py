#!/usr/bin/python3

import os, sys, re
from Bio import SeqIO

###CONFIGS###
# path_to_misas = "ANN_misas/edit_misas/"
# path_to_contigs = "ANN_contigs/"

path_to_misas = '/net/nfs/condor.sci.unh.edu/raid2/hcgs/joseph7e/anguina_genome/Sample_AT/N50_contig_files'
path_to_contigs = '/net/nfs/condor.sci.unh.edu/raid2/hcgs/joseph7e/anguina_genome/Sample_AT/N50_contig_files'

#path_to_sams = "NEM_sams/"
# outputfile = open("matching_ssr", 'w')
# outputfile.writelines("SAMPLE\tID\tSSR_type\tSSR\tSIZE\tSTART\tEND\tUPSTREAM\tDOWNSTREAM\tVERBOSE_ID\n")
flank = 30
how_many_common = 3
output_dir = "microstats/"

#IDEALY ADD A FINAL COLUMN WITH HETEROZYGOUS IDENTITY

def extract_file_names(search_term, path_to_dir):
    '''extracts a set of file names that contain a search term and returns it as a list'''
    filenames = []
    file_list = os.listdir(path_to_dir)
    for file_name in file_list:
        if search_term in file_name:
           filenames.append(path_to_dir + '/' + file_name)
    return filenames

misa_list = extract_file_names(".misa", path_to_misas)
contig_list = extract_file_names(".fasta", path_to_contigs)
#sam_list = extract_file_names(".sam", path_to_sams)

#print (misa_list, contig_list, sam_list)
print (misa_list)

contigs_misa = {}
for contig in contig_list:
    finder_tmp = contig.split("/"); finder, tmp = finder_tmp[-1].split(".")
    for misa in misa_list:
        if finder in misa and finder in contig:
            contigs_misa[contig] = misa
# for sam in sam_list:
#     tmp, finder_tmp = sam.split("/"); finder, tmp = finder_tmp.split(".")
#     for contig in contig_list:
#         for misa in misa_list:
#             if finder in misa and finder in contig and finder in sam:
#                 contigs_misaNsam[contig] = [misa,sam]
print ("\tCONTIG_FASTAS\t\t\tMICROSAT_FILES")
for c, m in contigs_misa.items():
    print (c, "--->", m)




def misa_extractor(line):
    sample, id, ssr_nr, ssr_type, ssr, size, start, end = line.split("\t")
    return (id, ssr_type, ssr, size, int(start), int(end))

def rev_comp(seq):
    """Reverses, complements and returns sequence"""
    rev_seq = seq[::-1]
    compliment_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    rev_comp_seq = ''
    for nuc in rev_seq:
        if nuc in ['A', 'T', 'G', 'C']:
            rev_comp_seq += compliment_dict[nuc]
    return rev_comp_seq

def extract_region(fasta, header_search, start, stop):
    """ PURPOSE: Extracts a specified region from a fasta using header_search, coordinates, and direction(+-),
                        and writes to append able output
        USAGE: extract_region(fasta, header_search, start, stop, output_name, direction, header_addition)"""
    for seq_record in SeqIO.parse(fasta, "fasta"):
        if header_search in str(seq_record.id):
            my_seq = (str(seq_record.seq[start-1:stop]))
    begin_flank = my_seq[:30]; end_flank = my_seq[-30:]; ssr = my_seq[30:-30]
    return (begin_flank, end_flank, ssr)

for c, m in contigs_misa.items():
    finder_tmp = c.split("/"); SAMPLE, tmp = finder_tmp[-1].split(".")
    outputfile = open(output_dir + SAMPLE + ".joe", 'w')
    outputfile.writelines("SAMPLE\tID\tSSR_type\tSSR\tSIZE\tSTART\tEND\tUPSTREAM\tDOWNSTREAM\tVERBOSE_ID\n")
    with open(m, 'r') as m:
        for line in m:
            if "SSR nr." in line:
                continue
            else:
                CONTIG, SSR_type, SSR, SIZE, START, END = misa_extractor(line)
                start = START - flank; end = END + flank
                UPSTREAM, DOWNSTREAM, VERBOSE_SSR = extract_region(c,CONTIG,start,end)
                print(SAMPLE+"\t"+CONTIG+"\t"+SSR_type+"\t"+SSR+"\t"+SIZE+"\t"+str(START)+"\t"+str(END)+"\t"+UPSTREAM+"\t"+DOWNSTREAM+"\t"+VERBOSE_SSR+"\n")
                outputfile.writelines(SAMPLE+"\t"+CONTIG+"\t"+SSR_type+"\t"+SSR+"\t"+SIZE+"\t"+str(START)+"\t"+str(END)+"\t"+UPSTREAM+"\t"+DOWNSTREAM+"\t"+VERBOSE_SSR+"\n")

def do_blast(command_list):
    '''performs blasts and returns the results'''
    sp = subprocess.Popen(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    blast = sp.communicate()
    blast = str(blast)
    return (blast)
