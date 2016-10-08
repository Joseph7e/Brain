#!/usr/bin/python3

import sys, os

input_dir = sys.argv[1]
output_dir = sys.argv[2]

replacer = ">0"


def extract_file_names(search_term, path_to_dir):
    '''extracts a set of file names that contain a search term and returns it as a list'''
    filenames = []
    file_list = os.listdir(path_to_dir)
    for file_name in file_list:
        if search_term in file_name:
           filenames.append(file_name)
    return filenames


fasta_files = extract_file_names(".fasta",input_dir)

for file in fasta_files:
    output = open(output_dir + file, 'w')
    count = 0
    with open(file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                line = line.replace(replacer, ">" + str(count))
                count += 1
            output.writelines(line)
print ("ALL DONE")
