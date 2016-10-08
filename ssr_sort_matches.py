#!/usr/bin/python3


import sys, os


file_list = os.listdir()

input_file = sys.argv[1]
output_file = open("sorted_polymorphic", 'w')
output_file2 = open("sorted_length", 'w')

block = ''
block_list = []
flag = False

with open(input_file, 'r') as f:
    for line in f:
        if line.startswith("#"):
            if block_list:
                if flag:
                    output_file2.writelines(block + "#\n")
                block_set = list(set(block_list))
                #print (block_list, '###',block_set, '######')
                if len(block_set) > 2:
                    output_file.writelines(block + "#\n")
            block = ''
            block2 = ''
            block_list = []
            flag = False
        else:
            block += line
            line_number, sample, node, ssr_type, ssr, length, start, stop, front, back, ssr_verbose, front_r, back_r, ssr_r_verbose = line.split("\t")
            block_list.append(length)
            if ssr_type != "c":
                tmp, total_repeats = ssr.split(")")
                if int(total_repeats) >= 8:
                    flag = True
            else:
                temper = ssr.split(")")
                if int(temper[-1]) >= 8:
                    flag = True

# header_ids = {'A02_N50':0, 'A12_N50':0, 'A3_N50':0, 'A4_N50':0, 'A6_N50':0, 'A7_N50':0, 'B02_N50':0, 'B03_N50':0, 'B1_N50':0, 'B7_N50':0, 'E9_N50':0, 'G2_N50':0, 'H2_N50':0}
# atlantic = ['A02_N50', 'A12_N50', 'A3_N50', 'A6_N50', 'A7_N50', 'B02_N50', 'B03_N50', 'B1_N50']
# pacific = ['A4_N0', 'B7_N50', 'E9_N50', 'G2_N50', 'H2_N50']








