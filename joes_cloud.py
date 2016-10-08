#!/usr/bin/python3

import sys
import mycloud

convert = True

string_match = ''

string_unmatch = ''

if convert:
    with open(sys.argv[1], 'r') as h:
        for line in h:
            line = line.rstrip()
            id, m, u = line.split('\t'); id = id + " "
            new_id_m = (id * int(m))
            new_id_u = (id * int(u))
            if m > u:
                string_match += new_id_m
            if u > m:
                string_unmatch += new_id_u
else:
    str = ''
    file = open(sys.argv[1], 'r')

    for line in file:
        str += line.rstrip()# + ","
    mycloud.generate_cloud(str,'word_out')

mycloud.generate_cloud(string_match, sys.argv[1] + "_enriched_matched_cloud")
mycloud.generate_cloud(string_unmatch, sys.argv[1] + "_enriched_unmatched_cloud")