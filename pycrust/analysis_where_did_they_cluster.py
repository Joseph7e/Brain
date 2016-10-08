#!/usr/bin/python3
#author: Joseph Sevigny

import re, sys
#input = '/home/genome/joseph7e/gene_16s_redo_all/complete_16S_uclust/analysis_what_matched_where/just_headers_aligned_16S_best_95_uc.fasta'
#input= '/home/genome/joseph7e/gene_16s_redo_all/complete_16S_uclust/analysis_what_matched_where/just_headers_aligned_amplified_16S_95_uc.fasta'

input = sys.argv[1]
input_uclust = open(input, 'r')



#output = open('what_matched_where_amplified_95.tsv', 'w')
output = open(input+'_matching.tsv', 'w')

gcf_dictionary = {}

for line in input_uclust.readlines():
    GCF = re.findall(r"(GCF_[0-9]*)",line)[0]
    cluster_number = re.findall(r">([0-9]*)|",line)[0]
    if GCF in gcf_dictionary.keys():
        if cluster_number in gcf_dictionary[GCF]:
            continue
        else:
            gcf_dictionary[GCF].append(cluster_number)
    else:
        gcf_dictionary[GCF] = []
        gcf_dictionary[GCF].append(cluster_number)
        #print (list(cluster_number))

single_count = 0
dub_count = 0
greater_count = 0

for key, value in gcf_dictionary.items():
    output_list = []
    output_list.append(key)
    if len(value) > 1:
        value = sorted(value,key=int)
        for thing in value:
            output_list.append(thing)
        if len(value) == 2:
            dub_count+=1
        if len(value) > 2:
            greater_count += 1
    else:
        output_list.append(value[0])
        single_count+=1
    output.writelines("\t".join(output_list))
    output.writelines('\n')
print (input)
print (single_count,dub_count,greater_count)
print (single_count+dub_count+greater_count)