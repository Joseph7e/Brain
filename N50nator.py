#!/usr/bin/python3

import sys, re


contig_file = sys.argv[1]
read_contig_file = open(contig_file, 'r')

N50_contig_file_name = sys.argv[2]
write_N = open(N50_contig_file_name, 'w')

cutoff = sys.argv[3]

cutoff = int(cutoff)

for line in read_contig_file.readlines():
    if ">" in line:
        node_length = re.findall(r"th_(\w*)_c", line)
        length = 0
        for i in node_length:
            length += int(i)
            if length < cutoff:
                quit()
            else:
                continue
    write_N.writelines(line)

print ("all_done")