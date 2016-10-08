#!/usr/bin/python3

input = open("final_download_list",'r')

output = open("amino_download_list",'w')

for line in input.readlines():
    if "GCF" in line and '.gff.gz' not in line:
        line = line.replace(".fna.gz", ".faa.gz")
        output.writelines(line)


