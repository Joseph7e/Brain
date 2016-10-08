#!/usr/bin/python3

#Author: Joseph Sevigny
#Affiliation: HCGS GEN 812
#Purpose: Run Busco, determine sequence overlap, extract sequences, and align for specified group.
#USAGE: ./busco_dealer.py <sample_commonality>

#nohup python3 ~/BUSCO/BUSCO_v1.1b1.py -l ~/BUSCO_v1.1b1/eukaryota/ -m genome -o Sample_NEM-H2 -in ~/smithsonian_reads/smithsonian/Project_Francesca_NemeAnnaNema_RP/Sample_NEM-H2/both_spades_output_error_correction/contigs.fasta &

import sys, subprocess, os, re
from Bio import SeqIO


            ####CONFIGS####
sample_name_commonality = sys.argv[1]
cutoff_value = int(sys.argv[2])

BUSCO = False
output_dir = "busco_dealer_common_SCP/"
nuc_output_dir = output_dir + "nucleotides/"
amino_output_dir = output_dir + "amino_acids/"
alignments_dir = output_dir + "alignments/"
sample_dir = output_dir + sample_name_commonality

if not os.path.exists(output_dir): #Create an output path for all created fastas and alignments
    os.mkdir(output_dir)
if not os.path.exists(nuc_output_dir):
    os.mkdir(nuc_output_dir)
if not os.path.exists(amino_output_dir):
    os.mkdir(amino_output_dir)
if not os.path.exists(alignments_dir):
    os.mkdir(alignments_dir)

                ###ALL THE DEFINITIONS USED THROUGHOUT####

def extract_file_names(search_term):
    '''extracts a set of file names that contain the
    specified file names and returns it as a list
    '''
    filenames = []
    file_list = os.listdir("./")
    for file_name in file_list:
        if search_term in file_name:
            filenames.append(file_name)
    return filenames

def gb2fasta(gb, fasta_name):
    input = open(gb, "rU")
    output = open(fasta_name, "w")
    sequences = SeqIO.parse(input, "genbank")
    SeqIO.write(sequences, output, "fasta")
    input.close()
    output.close()
    #print ("Converted {} records".format(count))

def get_gene_location(gb):
    gene_start = 0; gene_stop = 0
    for line in open(gb, 'r'):
        if "join" in line:
            continue
        elif "CDS" in line:
            line = line.replace("..", ":")
            if "complement" in line:
                line = line.replace(")", ":"); line = line.replace("(", ":")
                tmp1, start, stop, tmp2 = line.split(":")
                gene_start = int(start); gene_stop = int(stop)
            else:
                line = line.replace(" ", ""); line = line.replace("S", ":")
                tmp1, start, stop = line.split(":")
                gene_start = int(start); gene_stop = int(stop)
    return gene_start, gene_stop

def extract_region(fasta, start, stop, output_name, header_addition):
    with open(output_name, "a") as my_fasta:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            my_fasta.write(">" + header_addition + str(seq_record.id) + "\n")
            my_fasta.write(str(seq_record.seq[start-26:stop+24] + "\n"))

            ####RUN BUSCO IF YOU WOULD LIKE!!!####
if BUSCO:
    print("Running busco on all files, this is going to take a long time . . . ")
    path_to_dir = "~/smithsonian_reads/smithsonian/Project_Francesca_NemeAnnaNema_RP/"
    path_to_contigs = "/both_spades_output_error_correction/contigs.fasta"
    busco_output_names = ["Sample_NEM-A4", "Sample_NEM-A6"] #Retrieve list of all sample names from smithsonian folder
    def do_busco(command_list):
       sp = subprocess.Popen(command_list)
       blast = sp.communicate()
       #blast = blast.replace(r"\t", ":")
    for sample in busco_output_names:
       command = ['python3', '~/BUSCO/BUSCO_v1.1b1.py', '-l', '~/BUSCO/eukaryota/', '-m', 'genome', '-o', sample, '-in', path_to_dir + sample + path_to_contigs]
       print (command)
       do_busco(command)
    sys.exit()

                ####SEARCH FOR COMPLETE SCGs####

print ("Parsing the busco directories {}* for complete single copy orthologs ".format(sample_name_commonality))

busco_run_filenames = extract_file_names(sample_name_commonality) #grab all the busco runs

#Create a dictionary of scp_tags(keys) and corresponding list of sample names(values)
SCP_dictionary = {}
for current_sample_run in busco_run_filenames:
    current_scp_list = os.listdir(current_sample_run + "/selected/") #grabs the list of best complete scp
    temp_list = []
    for scp_name in current_scp_list:
        scp_name = scp_name[:-4]
        temp_list.append(scp_name)
    for scp in temp_list:
        if scp in SCP_dictionary:
            SCP_dictionary[scp].append(current_sample_run)
        else:
            SCP_dictionary[scp] = []
            SCP_dictionary[scp].append(current_sample_run)

#Create a new dictionary containing only the common SCPs
only_good_ones = {}
for scp, sample_name in SCP_dictionary.items():
    if len(sample_name) >= cutoff_value:
        only_good_ones[scp] = sample_name

print ("The following SCPs are the most common in your selected directories", only_good_ones)


#Enter directory for the sample and enter selected/ --> retrieve list of all complete genes - extension (.out)

#Enter into augustsus_proteins/ using list of complete genes extract the amino acid sequence

                ####CREATE AMINO ACID and NUCLEOTIDE MULTIFASTAS####

for key, value in only_good_ones.items():
    write_sequence_file = open(amino_output_dir + key + "_amino_sequences", 'w')
    for run in value:
        giff = open(run + "/augustus_proteins/" + key + ".fas", 'r')
        flag = 0
        for line in giff.readlines():
            line = line.replace('sequence', '')
            if line[0] == ">" and "g1" in line:
                new_head = ">" + run + line.replace(">", "_")
                write_sequence_file.writelines(new_head)
                flag = 1
            elif line[0] == ">" and "g2" in line:
                new_head = line.rstrip() + "_" + run + "\n"
                write_sequence_file.writelines(new_head)
                flag = 2
            elif line[0] == ">" and "g3" or "g4" in line:
                new_head = line.rstrip() + "_" + run + "\n"
                write_sequence_file.writelines(new_head)
                flag = 3
            if flag == 1 and line[0] != ">":
                write_sequence_file.writelines(line)
            if flag == 2 and line[0] != ">":
                write_sequence_file.writelines(line)
            if flag == 3 and line [0] != ">":
                write_sequence_file.writelines(line)


#use the list with matching hits to extract nucleotide sequences from genbank file

#convert all genbank files into fastas, determine start and stop location of scp, and extract from file
for key, value in only_good_ones.items():
    for run in value:
        gb_file = run + "/gb/" + key + ".raw.gb"
        fasta_sequence = nuc_output_dir + key + "_full_nuc_sequences.fasta"
        multi_fasta = nuc_output_dir + key + "_nucleotide.fasta"
        gb2fasta(gb_file, fasta_sequence)
        current_fasta = open(fasta_sequence, "r")
        gene_start, gene_stop = get_gene_location(gb_file)
        print(gene_start, gene_stop)
        extract_region(current_fasta, gene_start, gene_stop, multi_fasta, run)

#when they are all complete, create an alignmnet with muscle to comapre the genes

def muscle_alignment(input_multi_fastas, output_name):
    '''
    example_command = ["muscle", "-clw", "-in", "amino_fasta", "-out", "output_amino"]
    '''
    command = ["muscle", "-clw", "-in", input_multi_fastas, "-out", output_name]
    sp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    alignment = sp.communicate()


for key, value in only_good_ones.items():
    # amino_alignment_file_list = []
    # nucleotide_alignment_file_list = []
    multi_fasta = nuc_output_dir + key + "_nucleotide.fasta"
    amino_fasta = amino_output_dir + key + "_amino_sequences"
    output = alignments_dir + key +  "_alignment.clw"
    output_amino = alignments_dir + key + "_protein_alignment.clw"
    muscle_alignment(multi_fasta, output)
    muscle_alignment(amino_fasta, output_amino)

