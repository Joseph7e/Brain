#!/usr/bin/python3

import sys, re
contig_file = sys.argv[1]
blob_cov = sys.argv[2]

final_dict = {}

output = open('node_coverage_gc.txt', 'w')

for line in open(blob_cov):
    if line.startswith("#"):
        continue
    else:
        header, tmp, coverage = line.split("\t")
        nodes = re.findall("NODE_(\w*)_l", header)
        node = nodes[0]
        coverage = coverage.rstrip()
        final_dict[node] = [coverage]

def read_fasta(file_name):
    """reads a fasta file and returns a dictionary with headers as keys and sequences as values"""
    cfo = open(file_name, 'r')
    fasta_dict = {}
    for line in cfo.readlines():
        if line[0] == ">":
            header = line
        else:
            if header not in fasta_dict:
                fasta_dict[header] = line.rstrip()
            else:
                fasta_dict[header] += line.rstrip()
    return fasta_dict

print("Creating a dictionary of contig headers and sequences")
contig_dict = read_fasta(contig_file) #Uses the read_fasta definition to create a LARGE contig dictionary file

GC_count_list = []
Coverage_list = []
Length_list = []
#output = open(node_details_dir + sample_name + "_" + str(length_cutoff) + "_contig_details", 'w')
#output.writelines("#USE --> sort -t $'TAB' -nr -k3 details.file > sorted.file , To sort the file by the specified column (replace '3' with '2' or '4' for others)\n#\n#NODE_NAME\tCOV\tLength\tGC_Percentage\n")

#print("Determining node details, saving to output", node_details_dir)
for head, seq in sorted(contig_dict.items()): #sorts the contig nodes alphabetically
    #node_length = re.findall(r"th_(\w*)_c", head)
    #node_length = ''.join(node_length)
    #Length_list.append(int(node_length))
    length = len(seq)
    node = re.findall(">NODE_(\w*)_l", head) #Determines node name
    node = ''.join(node)
    #coverage = re.findall("_cov_(.*)_ID", head) #pulls coverage from the header, assumes spades output
    #coverage = ''.join(coverage);coverage = float(coverage)
    G = seq.count('G'); C = seq.count('C') #creates a count for each nucleotide in each of the expanded sequences
    GC_content = ((G+C)/len(seq)*100)
    GC_count_list.append(GC_content)
    #Coverage_list.append(coverage)
    #output.writelines("{0}\t{1:.2f}\t{2}\t{3:.2f}\n".format(node,coverage,node_length,GC_content))
    final_dict[node].append(GC_content)
    final_dict[node].append(length)

for key, value in final_dict.items():
    output.writelines(key+'\t'+str(value[2])+'\t'+str(value[0])+'\t'+str(value[1])+'\n')