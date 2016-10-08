#!/usr/bin/python3

#Purpose: Take a reference fasta with variant vcf and prodice a new fasta file. 
		#Take this file and gff and create a genbank file.

import sys

#arguments
reference_fasta = sys.argv[1]
variant_file = sys.argv[2]
#gff_file = sys.argv[3]


contig_order = []

def create_contig_dict(filename):
    """
    Quick and dirty creation of contigs dictionary from file.
    """
    contigs = {}
    with open(filename, "r") as f:
        header = ""
        seq = ""
        for line in f:
            if line[0] == ">":
                if header:
                    contigs[header] = seq
                header = line[1:].rstrip()
                contig_order.append(header.replace('>',''))
                seq = ""
            else:
                seq += line.rstrip()
        contigs[header] = seq
    return contigs

contig_dict = create_contig_dict(reference_fasta)


addition = 0
current_contig = ''
previous_range_start = 0
previous_range_end = 0


for line in open(variant_file, "r").readlines():
    if line[0] == '#':
        continue
    elements = line.rstrip().split("\t")
    name = elements[0] # name of the contig / chromosome SNP found on
    pos = elements[1] # position of SNP in contig
    ref = elements[3] # reference basepair(s)
    alt = elements[4] # SNP basepair(s) at same location
    if ',' in alt:
        alt = alt.split(',')[0]
    qual = float(elements[5]) # SNP quality, might want to filter this
    if current_contig != name: # make sure to reset counts when hitting a new contig
        print ('Working on contig', [name],'\n\n\n')
        current_contig = name
        addition = 0; previous_range_start = 0; previous_range_end = 0

    ####set reference start and stops
    contig_length = len(contig_dict[name])
    start_pos = addition + int(pos) -1; end_pos = start_pos + len(ref)
    new_end = start_pos + len(alt) # troubleshooting

    ### Avoid multiple start from non-existent reference location
    if int(pos) > previous_range_start and int(pos) < previous_range_end:
        print ("SNP OUT OF RANGE\n\n\n\n\n")
        continue
    previous_range_start = int(pos); previous_range_end = int(pos)+len(ref)

    print ("original snp and flanks")
    print (contig_dict[name][start_pos-5:start_pos], contig_dict[name][start_pos:end_pos], contig_dict[name][end_pos:end_pos+5])
    print ("contig,pos,ref,alt,qual")
    print (name,pos,ref,alt,qual)
    if contig_dict[name][start_pos:end_pos] != ref:
        print ('SNP identification and and reference location do not match')
        new = contig_dict[name][:start_pos] + alt + contig_dict[name][end_pos:]
        contig_dict[name] = new
        print (start_pos, end_pos)
        print (contig_dict[name][start_pos-5:start_pos], contig_dict[name][start_pos:new_end], contig_dict[name][new_end:new_end+5])
        addition += (len(alt) - len(ref))#track correct start locations for shifting reference
        sys.exit()
    #print (contig_dict[name][start_pos:end_pos])


    new = contig_dict[name][:start_pos] + alt + contig_dict[name][end_pos:]
    contig_dict[name] = new
    print (contig_dict[name][start_pos-5:start_pos], contig_dict[name][start_pos:new_end], contig_dict[name][new_end:new_end+5])
    print (start_pos, end_pos)
    addition += (len(alt) - len(ref)) #track correct start locations for shifting reference
    print (addition)
    print ()
#print (contig_dict['chrXVI'][946510:946511])

outfile = open("it_worked.fasta","w")


for c in contig_order:
    outfile.writelines(">"+c+"\n"+contig_dict[c]+"\n")