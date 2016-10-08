#!/usr.bin/python3

#Author: Joseph7e
#Purpose: Parse a vcf for haplotype information (used w/ various ploidy levels)

import sys, re
vcf = sys.argv[1]

#outname = open(vcf+'parsed','w')

all_ratios = []
all_second_alts = []
max_alts = 1; alt_count = 0
total_count = 0
high_ratio = 0; low_ratio = 9999999
low_qual_cutoff = 50

contig_dict = {}


for line in open(vcf,'r').readlines():
    if line.startswith('#'):
        continue
    else:
        elements = line.split('\t')
        contig = elements[0]
        location = elements[1]
        ref = elements[3]
        alt = elements[4]
        alts = []
        if ',' in alt:
            alts = alt.split(',')
        qual = elements[5]
        filter = elements[6]
        info = elements[7]
        better_info = elements[9]
        ratios = better_info.split(':')[1].split(',')
        if filter == "LowQual" or float(qual) < low_qual_cutoff or int(ratios[0]) < 3 or int(ratios[1]) < 3:
            continue
        else:
            total_count += 1
            if alts:
                alt_count +=1
                if max_alts >= len(alts):
                    continue
                else:
                    max_alts = len(alts)
                count = 1
                alt_total = 0
                alt = ''
                for alter in alts:
                    alt += alter+','
                    alt_total += int(ratios[count])
                    count +=1
                ratio = (alt_total/(alt_total+int(ratios[0]))*100)
                alt = alt[:-1]
                ratio_tmp = [ratios[0],alt_total]
            else:
                ratio_tmp = better_info.split(':')[1].split(',')
            if '0' in ratio_tmp:
                ratio = 0
            if int(ratio_tmp[1]) > int(ratio_tmp[0]) and '0' not in ratio_tmp:
                ratio = int(ratio_tmp[0])/(int(ratio_tmp[1])+int(ratio_tmp[0]))
            else:
                ratio = int(ratio_tmp[1])/(int(ratio_tmp[0])+int(ratio_tmp[1]))
            length = re.findall(r'length_([0-9]*)_',contig)[0]
            contig = re.findall(r'(NODE_[0-9]*)_',contig)[0]
            if ratio > high_ratio:
                high_ratio = ratio
            if ratio < low_ratio:
                low_ratio = ratio
            print ('\nContig ID = {}\nLocation(quality) = {}({}) , Reference(depth) = {}({}), Alternate(depth) = {}({})\nRatio = {}'.format(contig, location, qual, ref, str(ratio_tmp[0]),alt,str(ratio_tmp[1]), str(int(ratio*100))+'%'))
            all_ratios.append(ratio)
            if contig in contig_dict.keys():
                contig_dict[contig][0] += 1
            else:
                contig_dict[contig] = [1, length]
print ("\n\nTotal number of high quality snps = ", total_count)
print ('Average ratio of alternate allele depths = {0:.2f}%'.format((sum(all_ratios)/len(all_ratios))*100))
print ('Maximum number of alternate alleles = ', max_alts, '\nTotal number of alleles with more than one alt = {}  ({:.2f}% of total snps)'.format(alt_count,(alt_count/total_count)*100))
print ('High ratio = {:.2f}%\tLow ratio = {:.2f}%'.format(high_ratio*100,low_ratio*100))

for c, co in contig_dict.items():
    print (c + '\t' + str(co[0]) + '\t' + str(co[1]) + '\t' + str(round((co[0]/int(co[1])),4)))
