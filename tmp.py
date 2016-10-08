#!/usr/bin/python3

import os

cluster_dir = 'clusters_fixed/'
tbl_output = open('graph', 'w')
#cluster_files = os.listdir(cluster_dir)
cluster_files = ['data_100.tsv', 'data_99.tsv', 'data_98.tsv', 'data_97.tsv', 'data_96.tsv', 'data_95.tsv']


for line in open('all_round_comparsion_result.csv','r').readlines():
    flag = False
    file_name, total_expected, asigned, unassigned, ratio = line.split(',')
    ratio = ratio.rstrip()
    for cluster in cluster_files:
        if flag == False:
            with open(cluster_dir + cluster, 'r') as c:
                for name in c.readlines():
                    name = name.rstrip()
                    if name == file_name:
                        tbl_output.writelines(cluster[:-4] + '\t' + ratio + '\t' + total_expected + '\t' + name + '\n')
                        flag = True
        else:
            continue

