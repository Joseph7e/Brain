#!/usr/bin/python3

import os, subprocess

fastas_dir = "./"
output_alignments_dir = "alignments/"

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


new_fastas = extract_file_names("fa", fastas_dir)

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
    muscle_alignment(fastas_dir + current_fasta, output)
