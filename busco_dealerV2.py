#!/usr/bin/python3

#Author: Joseph Sevigny
#Affiliation: HCGS
#Purpose: Run Busco, determine sequence overlap, extract sequences, and align for specified group.
#USAGE: ./busco_dealer.py <sample_commonality> <Path to Contig Directory> <number in common(int)>

import sys, subprocess, os, re
from Bio import SeqIO

sample_name_commonality = sys.argv[1] #run_Sample_Nem
contig_dir = sys.argv[2] #~/smithsonian/smithsonian/ProjectFrancesca_RP/
cutoff_value = int(sys.argv[3]) #13

##note: may need to adjust a few downstream pipelines that grab a range of the file_name if normal busco convention is not used


        ######CONFIGS#######

path_to_busco_runs = "./"
path_to_to_contigs = "/both_spades_output_error_correction/contigs.fasta" #path from contig_dir to contig.fasta
RUN_BUSCO = "NO" #set to "NO" to bypass running busco for each sample or "YES" to run busco

output_dir = "b_dealer_V3_/"
nuc_output_dir = output_dir +  sample_name_commonality + "nucleotides/"
amino_output_dir = output_dir + sample_name_commonality + "amino_acids/"
alignments_dir = nuc_output_dir +  "alignments/"

if not os.path.exists(output_dir): #creates a directory hierarchy for files
    os.mkdir(output_dir)
if not os.path.exists(nuc_output_dir):
    os.mkdir(nuc_output_dir)
if not os.path.exists(amino_output_dir):
    os.mkdir(amino_output_dir)
if not os.path.exists(alignments_dir):
    os.mkdir(alignments_dir)


        ######RUN BUSCO########

def run_busco(contigs):
    print ("ADD ME")
    #nohup python3 ~/BUSCO/BUSCO_v1.1b1.py -l ~/BUSCO_v1.1b1/eukaryota/ -m genome -o Sample_NEM-H2 -in ~/smithsonian_reads/smithsonian/Project_Francesca_NemeAnnaNema_RP/Sample_NEM-H2/both_spades_output_error_correction/contigs.fasta &


if RUN_BUSCO == "YES":
    print ("Running busco")

        ######DETERMINE ALL COMMON and COMPLETE SINGLE COPY GENES#####

print ("\nParsing the busco directories {}* for complete single copy orthologs ".format(sample_name_commonality))


def extract_file_names(search_term, path_to_dir):
    '''extracts a set of file names that contain a search term and returns it as a list'''
    busco_run_filenames = []
    file_list = os.listdir(path_to_dir)
    for file_name in file_list:
        if search_term in file_name:
            busco_run_filenames.append(file_name)
    return busco_run_filenames


busco_run_filenames = extract_file_names(sample_name_commonality, path_to_busco_runs) #Creates a list of all the busco runs

contig_file_names = [] #creates a matching list of contig file names
for run in busco_run_filenames:
    sample = run[4:-3]
    contig_loci = extract_file_names(sample, contig_dir)
    contig_file_name = contig_loci[0] + path_to_to_contigs
    contig_file_names.append(contig_dir + contig_file_name)


#Create a dictionary of scp_tags(keys) and corresponding list of sample names(values)

def extract_complete_busco_names(sample):
    '''Parses buscos full_table of all scp genes and returns only the identity of the complete ones as a list'''
    sample_name = sample[4:]
    complete_genes = []
    try:
        for line in open(sample + "/full_table_" + sample_name):
            if "Complete" in line:
                identity, complete, node, tmp1, tmp2, tmp3, tmp4 = line.split("\t")
                complete_genes.append(identity)
    except FileNotFoundError:
        sample_name = sample[11:]
        for line in open(sample + "/full_table_" + sample_name):
            if "Complete" in line:
                identity, complete, node, tmp1, tmp2, tmp3, tmp4 = line.split("\t")
                complete_genes.append(identity)
    return complete_genes




SCP_dictionary = {}
for current_sample_run in busco_run_filenames:
    current_scp_list = os.listdir(current_sample_run + "/augustus/") #grabs the list of all scps
    temp_list = []
    completes = extract_complete_busco_names(current_sample_run)
    for scp_name in current_scp_list:
        file = current_sample_run + "/augustus/" + scp_name
        try:
            scp_name, ext = scp_name.split(".")
        except ValueError:
            scp_name, ext, number = scp_name.split(".")
        if scp_name in completes:
            temp_list.append(scp_name)
        for scp in temp_list:
            if scp in file:
                if scp in SCP_dictionary:
                    SCP_dictionary[scp].append(file)
                else:
                    SCP_dictionary[scp] = []
                    SCP_dictionary[scp].append(file)

#Create a new dictionary containing only the common SCPs

only_good_ones = {}
for scp, sample_names in SCP_dictionary.items():
    sample_names = list(set(sample_names)) #make sure we only have unique sample_names
    keeping_track = []
    for sample_name in sample_names:
        if sample_name[11:18] in keeping_track:
            continue
        else:
            keeping_track.append(sample_name[11:18])
    if len(keeping_track) >= cutoff_value:
        only_good_ones[scp] = sample_names
print ("\nThe total number of genes common to at least {} of the samples is {} :\n".format(str(cutoff_value), len(only_good_ones)))
#####CREATE A MULTI-FASTA FOR EACH OF THE COMMON SCGS#####

print ("Extracting single_copy_genes from each sample and writing to multifasta output\nThis may take a bit ...\n")

def extract_region(fasta, header_search, start, stop, output_name, direction, header_addition):
    """ PURPOSE: Extracts a specified region from a fasta using header_search, coordinates, and direction(+-),
                        and writes to append able output
        USAGE: extract_region(fasta, header_search, start, stop, output_name, direction, header_addition)"""
    with open(output_name, "a") as my_fasta:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            if header_search in str(seq_record.id):
                my_fasta.write(">" + header_addition + str(seq_record.id) + "\n")
                my_seq = (str(seq_record.seq[start-1:stop]))
                if direction == "-":
                    rev_seq = rev_comp(my_seq)
                    my_fasta.write(rev_seq + "\n")
                else:
                    my_fasta.write(my_seq + "\n")

def rev_comp(seq):
    """Reverses, complements and returns sequence"""
    rev_seq = seq[::-1]
    compliment_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    rev_comp_seq = ''
    for nuc in rev_seq:
        if nuc in ['A', 'T', 'G', 'C']:
            rev_comp_seq += compliment_dict[nuc]
    return rev_comp_seq


for scp_name, sample_names in only_good_ones.items():
    print ("\n\n***** Currently working on single copy gene --> ", scp_name, "*****\n")
    count = 0
    pro_count = 0
    current_output = nuc_output_dir + scp_name + "_nucleotide_sequences.fasta"
    protein_output = amino_output_dir + scp_name + "_amino_acid_sequences.fasta"
    write_protein =  open(protein_output, 'w')
    for run in sample_names:
        tmpbegin, tmpscp = run.split("augustus/")
        scp_tmp_name, tmpext = tmpscp.split("out")
        protein_file = tmpbegin + "augustus_proteins/" + scp_tmp_name + "fas" + tmpext
        print (protein_file)
        print (run)
        with open(protein_file, 'r') as p:
            for seq_line in p:
                seq_line = seq_line.replace("sequence", "")
                seq_line = seq_line.replace(">", ">" + str(pro_count) + run[11:20] + "_" + scp_name)
                write_protein.writelines(seq_line)
                pro_count += 1
        #Could make below a function????
        with open(run, 'r') as f:
            for line in f.readlines(): #parse through the augustus files and extract start and stop regions
                if line[0] == "#":
                    continue
                else:
                    if "gene" in line:
                        header, tmp, gene, start, stop, tmp3, direction, tmp4, tmp5 = line.split("\t")
                        start = int(start); stop = int(stop)
                        if gene == "gene": #A double check to make sure we are capturing genes
                            for contig_file in contig_file_names:
                                sample, tmp, tmp2 = run.split("/")
                                if sample[4:-3] in contig_file:
                                    with open(contig_file) as c:
                                        print ("File = {}".format(run + "-->" + header))
                                        print ("Start = {}, Stop = {}  Length = {}, Direction = ({})".format(start, stop, stop-start, direction), "\n")
                                        extract_region(c, header, start, stop, current_output, direction, str(count) + run[11:20] + run[27:] +  "_" + str(start) + "-" + str(stop) + "(" + direction + ")")
                                        count += 1


#when they are all complete, create an alignment with muscle to compare the genes

def muscle_alignment(input_multi_fastas, output_name):
    '''
    example_command = ["muscle", "-clw", "-in", "amino_fasta", "-out", "output_amino"]
    '''
    command = ["muscle", "-clw", "-in", input_multi_fastas, "-out", output_name]
    sp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    alignment = sp.communicate()


new_files = extract_file_names(".fasta", nuc_output_dir)
for current_file in new_files:
    #amino_fasta = amino_output_dir + key + "_amino_sequences"
    output = alignments_dir + current_file +  "_alignment.clw"
    #output_amino = alignments_dir + key + "_protein_alignment.clw"
    muscle_alignment(current_file, output)
    #muscle_alignment(amino_fasta, output_amino)

#USE THE ALIGNMENT TO ADJUST SCP GENE FASTAS AND RERUN do_alignments.py
