#!/usr/bin/python

import sys, os, re

def extract_file_names(search_term, path_to_dir):
    '''extracts a set of file names that contain a search term and returns it as a list'''
    filenames = []
    file_list = os.listdir(path_to_dir)
    for file_name in file_list:
        if search_term in file_name:
           filenames.append(file_name)
    return filenames

file_list = extract_file_names(".txt", sys.argv[1])

output = open(sys.argv[2], 'w')

for file in file_list:
    genome = file.split('.')[0]
    count = 0
    with open(sys.argv[1] + file) as p:
        for line in p:
            if "......" in line:
                print (genome)
                try:
                    sequence = re.findall(r"\s([CATG]*)\n", line)[0]
                    if sequence == '':
                        print ("blank sequence")
                        continue
                    else:
                        count += 1
                        output.writelines(">" + str(count) + "_" + genome + '\n' + sequence + '\n')
                except IndexError:
                    print ("lost one -->", sequence)
                    continue
