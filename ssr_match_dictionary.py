#!/usr/bin/python3


import os, sys
match_list = os.listdir(sys.argv[1])
output_table = open(sys.argv[2], 'w')

for file in match_list:
    count = 0
    with open(sys.argv[1] + file, 'r') as f:
        for line in f:
            if "NEM" in line:
                count += 1
    file = file.replace("_matches", "")
    output_table.writelines(str(count) + "\t" + file + "\n")