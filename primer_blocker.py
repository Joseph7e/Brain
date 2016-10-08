#!/usr/bin/python3

import sys, subprocess




reference_fasta = '/home/genome/joseph7e/primer_block/reference_16S_amplified.fasta'
subject_sequence = 'anthrobacter_amplified_region.fasta'
kmer_size = 15

complete_sequence = ''

with open(subject_sequence,'r') as s:
    for line in s.readlines():
        if ">" not in line and 'A' in line:
            complete_sequence = line

kmer_dict = {}

count = 0

current_position_start = 0
print ('primer,','start :','exact,','exact-1,','exact-2,','exact-3,','exact-4,')

while current_position_start < (len(complete_sequence) - kmer_size):
    kmer = (complete_sequence[current_position_start:current_position_start+kmer_size])
    kmer_dict[current_position_start] = kmer
    current_position_start += 1
for start, k in kmer_dict.items():
    command = ['grep', '-c', k,reference_fasta]
    sp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    grep_score = str(sp.communicate()[0].decode('ascii'));grep_score = int(grep_score)
    if grep_score < 4:
        scores = []
        for i in range(6):
            if i != 0:
                command2 = ['agrep', '-'+str(i), '-c', k, reference_fasta]
                sp = subprocess.Popen(command2, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                score = str(sp.communicate()[0].decode('ascii')).replace('\n','')
                score = int(score)
                scores.append(score)
        print (k, start-1, ':', str(grep_score)+',', str(scores[0]) + ',', str(scores[1]) + ',', str(scores[2]) + ',', str(scores[3]))


