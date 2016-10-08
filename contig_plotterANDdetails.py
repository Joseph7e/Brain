#!/usr/bin/python3

#Author: Joseph Sevigny
#Affiliation: University of New Hampshire, Hubbard Center for Genome Studies
#Purpose: Create scatter plots to examine potential contamination in spades output contigs
#USAGE: ./contig_plotterANDdetails.py path/to/contigs.fasta Sample_name [contig_length_cutoff]

import sys, re, os, matplotlib.pyplot as pyplot

help = """
                                Welcome to contig_plotterANDdetails.py

        AUTHOR: Joseph Sevigny - Hubbard Center For Genome Studies
        Purpose: create scatter plots and table of contig details (examine potential contamination)
        Assumptions: standard SPADES assembler output with length and coverage in the headers
        USAGE: ./contig_plotterANDdetails.py <path/to/contigs.fasta> <Sample_name> [contig_length_cutoff]
        default cutoff value = 0

        Example Usage --> ./contig_plotterANDdetails.py files/ANN_B04_trimmed_contigs.fasta ANN-B04 1000

        """


#CONFIGS
try:
    contig_file = sys.argv[1] # path to contig file
    if contig_file == "-h":
        print (help)
        sys.exit()
except IndexError:
    print (help)
    sys.exit()
sample_name = sys.argv[2] #A name to label graphs and output files
length_cutoff = 0 # default length cutoff, can be changed with sys.argv[3]
try:
    length_cutoff = sys.argv[3]
except IndexError:
    pass

output_dir = "scatter_plots/"
if not os.path.exists(output_dir): #Creates an output dir for node details and scatter plots
    os.mkdir(output_dir)
node_details_dir = "scatter_plots/node_details/"
if not os.path.exists(node_details_dir):
    os.makedirs(node_details_dir)


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

#node detail lists
GC_count_list = []
Coverage_list = []
Length_list = []
output = open(node_details_dir + sample_name + "_" + str(length_cutoff) + "_contig_details", 'w')
output.writelines("#USE --> sort -t $'TAB' -nr -k3 details.file > sorted.file , To sort the file by the specified column (replace '3' with '2' or '4' for others)\n#\n#NODE_NAME\tCOV\tLength\tGC_Percentage\n")

print("Determining node details, saving to output", node_details_dir)
for head, seq in sorted(contig_dict.items()): #sorts the contig nodes alphabetically
    node_length = re.findall(r"th_(\w*)_c", head)
    node_length = ''.join(node_length)
    if int(node_length) >= int(length_cutoff): #filters only the nodes with a certain length
        Length_list.append(int(node_length))
        node = re.findall(">NODE_(\w*)_l", head) #Determines node name
        node = ''.join(node)
        coverage = re.findall("_cov_(.*)_ID", head) #pulls coverage from the header, assumes spades output
        coverage = ''.join(coverage);coverage = float(coverage)
        G = seq.count('G'); C = seq.count('C') #creates a count for each nucleotide in each of the expanded sequences
        GC_content = ((G+C)/len(seq)*100)
        GC_count_list.append(GC_content)
        Coverage_list.append(coverage)
        output.writelines("{0}\t{1:.2f}\t{2}\t{3:.2f}\n".format(node,coverage,node_length,GC_content))


def plotdata(x_list, y_list, title, subtitle, x_label, y_label, out_name, sub_output):
    """ takes a user defined set of data and creates a jpg graph"""
    #x_data = x_list
    #y_data = y_list
    figure = pyplot.figure()
    figure.suptitle(sample_name, fontsize=20, fontweight='bold')
    figure.subplots_adjust(top=0.85)
    ax = figure.add_subplot(111)
    ax.set_title(subtitle)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    pyplot.scatter(x_list,y_list, color="blue", marker=".") #may want to use marker="," if . is too big
    pyplot.savefig(output_dir + out_name + "_cutoff_" + sub_output + "_scatterplot.jpg")


print("Creating a scatter plot of GC content vs Coverage") #Read print for notes, running plotdata for all data
length_cutoff = str(length_cutoff)
plotdata(Coverage_list,GC_count_list, sample_name, "Contigs (>=" + length_cutoff +  ") : Coverage vs GC Percentage", "Coverage", "GC Content", sample_name, length_cutoff + "_GC")

print("Creating a scatter plot of Length vs GC content")
plotdata(Length_list,GC_count_list, sample_name, "Contigs (>=" + length_cutoff +  ") : Length vs GC Percentage", "Length", "GC Content", sample_name, length_cutoff + "_Len")

print("\nProcess complete, a total of {} contigs in the file greater than {}bp\nScatter plots and tables can be found in {}\nGoodbye".format(len(GC_count_list), length_cutoff, output_dir))
