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
                if len(block_list) == len(block_set):
                    output_file.writelines(block + "#\n")
            block = ''
            block2 = ''
            block_list = []
            flag = False
        else:
            block += line
            line_number, sample, node, ssr_type, ssr, length, start, stop, front, back, ssr_verbose, front_r, back_r, ssr_r_verbose = line.split("\t")
            block_list.append(length)
            if ssr_type != "c*":
                tmp, total_repeats = ssr.split(")")
                if int(total_repeats) >= 8:
                    flag = True
