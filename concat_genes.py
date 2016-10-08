#!/usr/bin/python3

import os, sys
from Bio import SeqIO
# from shutil import copyfile


#
# os.mkdir('housekeeping/' + sys.argv[1])
#
# for fasta in os.listdir(sys.argv[1]):
#     count = 0
#     with open(sys.argv[1] + fasta, 'r') as f:
#         for line in f:
#             if line.startswith(">"):
#                 count += 1
#
#         if count == 1:
#             copyfile(sys.argv[1] + fasta, 'housekeeping/' + sys.argv[1] + fasta)

#
output = open("full_table.tsv", 'w')

output2 = open("best_table.tsv", 'w')

dict = {} # genome: value

def list_dir_grab_id(dir):
    list = os.listdir(dir)
    for g in list:
        g = g[:13]
        if g in dict.keys():
            dict[g].append(dir[13:])
        else:
            dict[g] = [dir[13:]]



list_dir_grab_id('gene_extract_initiation_factor')
list_dir_grab_id('gene_extract_gyrase')
list_dir_grab_id('gene_extract_ctp_synthase')
list_dir_grab_id('gene_extract_ribosomal_protein_S12')
list_dir_grab_id('gene_extract_ribosomal_protein_S15')


for key, value in dict.items():
    output.writelines(key + "\t" +  str(value) + '\n')
    if len(value) == 5:
        output2.writelines(key + '\n')





#gene_count = 13
#
# if_2 = os.listdir()
# gyrase = os.list()
# ctp =
# s12 =
# s15 =
#
#
#
# def extract_file_names(search_term, path_to_dir):
#     '''extracts a set of file names that contain a search term and returns it as a list'''
#     filenames = []
#     file_list = os.listdir(path_to_dir)
#     for file_name in file_list:
#         if search_term in file_name:
#            filenames.append(file_name)
#     return filenames
#
# file_names = extract_file_names(".fa", "./")
#
# output = open("concat_sequences.fa", 'w')
#
# print (file_names)
#
# headers = []
#
# for seq_record in SeqIO.parse(file_names[5], "fasta"):
#     print (seq_record.id[5:15])
#     headers.append(seq_record.id[5:30])
#
#
# for header in headers:
#     count = 0
#     output.writelines(">" + header + "\n")
#     for file in file_names:
#          for seq_record in SeqIO.parse(file, "fasta"):
#              if header in seq_record.id:
#                  output.writelines(str(seq_record.seq))
#                  count += 1
#     output.writelines("\n")
#     if count < gene_count:
#         print (header + "\t<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
#     if count > gene_count:
#         print (header + "\t>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
#     print (count)
#     # for seq_record in SeqIO.parse(contig, "fasta"):
#     #                 if current_node in str(seq_record.id):
#     #                     cfo.write(">" + contig + str(seq_record.id) + "\n")
#     #                     cfo.write(str(seq_record.seq) + "\n")