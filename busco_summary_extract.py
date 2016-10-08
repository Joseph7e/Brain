#!/usr/bin/python3

import os

path_to_summary = "/short_summary_Sample"
path_to_summary2 = "/short_summary"

output = open("busco_summary_for_all.xml", 'w')

line0 = 'Sample\t'
line1 = 'Complete_Single-Copy_BUSCOs\t'
line2 = 'Complete_Duplicated_BUSCOs\t'
line3 = 'Fragmented_BUSCOs\t'
line4 = 'Missing_BUSCOs\t'
line5 = 'Total_BUSCO_groups_searched\t'

def extract_file_names(search_term, path_to_dir):
    '''extracts a set of file names that contain a search term and returns it as a list'''
    busco_run_filenames = []
    file_list = os.listdir(path_to_dir)
    for file_name in file_list:
        if search_term in file_name:
            busco_run_filenames.append(file_name)
    return busco_run_filenames

busco_samples = extract_file_names("Sample", "." )

buco_samples = busco_samples.sort()

for sample in busco_samples:
    line0 += sample + "\t"
    try:
        path = sample + path_to_summary + sample[10:]
        for line in open(path):
            if "Single-Copy" in line:
                tmp, number, title = line.split("\t")
                line1 += number + "\t"
            if "Duplicated" in line:
                tmp, number, title = line.split("\t")
                line2 += number + "\t"
            if "Fragmented" in line:
                tmp, number, title = line.split("\t")
                line3 += number + "\t"
            if "Missing" in line:
                tmp, number, title = line.split("\t")
                line4 += number + "\t"
            if "Total" in line:
                tmp, number, title = line.split("\t")
                line5 += number + "\t"
    except FileNotFoundError:
        path = sample + path_to_summary2 + sample[10:]
        for line in open(path):
            if "Single-Copy" in line:
                tmp, number, title = line.split("\t")
                line1 += number + "\t"
            if "Duplicated" in line:
                tmp, number, title = line.split("\t")
                line2 += number + "\t"
            if "Fragmented" in line:
                tmp, number, title = line.split("\t")
                line3 += number + "\t"
            if "Missing" in line:
                tmp, number, title = line.split("\t")
                line4 += number + "\t"
            if "Total" in line:
                tmp, number, title = line.split("\t")
                line5 += number + "\t"
output.writelines(line0 + "\n")
output.writelines(line1 + "\n")
output.writelines(line2 + "\n")
output.writelines(line3 + "\n")
output.writelines(line4 + "\n")
output.writelines(line5)

print ("done")