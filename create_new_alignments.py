#!/usr/bin/python3

import os, sys, subprocess
from Bio import SeqIO


busco_list = "../clusters.txt"

#busco_list will include something like this
# 38531   2,3,4,5,6,7,8,9,10,11,15
# 108683  22,10,18,25,5,16,13,4,21,9
# 35347   33,3,4,18,36,10,24,12,30,16,21
# 38329   1,2,3,5,7,9,10,11,12,13,14
# 53371   6,22,26,8,14,25,18,12,5,23
# 12273   1,2,11,23,45,9,20,16,4  41,18,37,47,31,39,24,32,29      25,22,29,26,43,15,7,8,21
# 15263   33,1,3,10,12,21,29,26,9,27,20
# 36637   34,23,26,22,41,3,16,36,33,31,44 46,24,39,15,7,32,28,11,21,30,9
# 46157   15,10,17,6,12,13,16,7,3,1,9
# 51307   13,43,44,16,5,30,34,40,1,12
# 101066  1,3,7,5,14,10,2,8,9,13
# 108369  9,6,20,12,8,13,19,4,2,1
# 46157   15,10,17,6,12,16,7,3,1,9
# 52229   19,5,14,27,2,3,25,12,1,20
# 55553   1,12,28,19,33,13,36     4,35,8,10,34,18,27,6,24,2,11
# 60479   18,23,32,25,16,17,24,29,38,19,4,31,39   20,33,2,15,36,10,13,3,27,22,1   18,23,32,25,16,17,24,29,38,19,4,31,39,20,33,2,15,36,10,13,3,27,22,1


fasta_dir = "original_files/"
search_term  = "fasta"

output_fastas_dir = "busco_new_fastas/"
output_alignments_dir = "busco_new_alignments/"

if not os.path.exists(output_fastas_dir):
    os.makedirs(output_fastas_dir)

if not os.path.exists(output_alignments_dir):
    os.makedirs(output_alignments_dir)


def extract_file_names(search_term, path_to_dir):
    '''extracts a set of file names that contain a search term and returns it as a list'''
    filenames = []
    file_list = os.listdir(path_to_dir)
    for file_name in file_list:
        if search_term in file_name:
            filenames.append(file_name)
    return filenames

old_fastas = extract_file_names(search_term, fasta_dir)

oldfastas_list_of_names = {}

for oldf in old_fastas:
    with open(busco_list, 'r') as bl:
        for line in bl:
            line = line.rstrip()
            try:
                busco_id, list = line.split("\t")
                list = list.split(",")
                if busco_id in oldf:
                    oldfastas_list_of_names[oldf + "1"] = list
            except ValueError:
                try:
                    busco_id, list1, list2 = line.split("\t")
                    list1 = list1.split(",")
                    list2 = list2.split(",")
                    if busco_id in oldf:
                        oldfastas_list_of_names[oldf + "1"] = [list1]
                        oldfastas_list_of_names[oldf + "2"] = [list2]
                except ValueError:
                    busco_id, list1, list2, list3 = line.split("\t")
                    list1 = list1.split(",")
                    list2 = list2.split(",")
                    list3 = list3.split(",")
                    if busco_id in oldf:
                        oldfastas_list_of_names[oldf + "1"] = [list1]
                        oldfastas_list_of_names[oldf + "2"] = [list2]
                        oldfastas_list_of_names[oldf + "3"] = [list3]

for oldfasta, list in oldfastas_list_of_names.items():
    with open(output_fastas_dir + oldfasta, 'w') as new_fasta:
        for seq_record in SeqIO.parse(fasta_dir + oldfasta[:-1], "fasta"):
            for digits in list:
                search = str(digits) + "NEM"
                if str(seq_record.id).startswith(search):
                    print (search)
                    print (seq_record.id)
                    new_fasta.write(">" + str(seq_record.id) + "\n")
                    new_fasta.write(str(seq_record.seq) + "\n")

new_fastas = extract_file_names("fasta", output_fastas_dir)

def muscle_alignment(input_multi_fastas, output_name):
    '''
    example_command = ["muscle", "-clw", "-in", "amino_fasta", "-out", "output_amino"]
    '''
    command = ["muscle", "-clw", "-in", input_multi_fastas, "-out", output_name]
    sp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    alignment = sp.communicate()

for current_fasta in new_fastas:
    output = output_alignments_dir + current_fasta + ".clw"
    print (current_fasta, output)
    muscle_alignment(output_fastas_dir + current_fasta, output)