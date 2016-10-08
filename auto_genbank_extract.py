#!/usr/bin/python3
import subprocess, sys, os

path_to_file = sys.argv[1]
count = 1


fasta_output_dir = "fastas/"
if not os.path.exists(fasta_output_dir): #creates a output path for all mitochondrial nodes
    os.mkdir(fasta_output_dir)
genbank_output_dir = "genbank/"
if not os.path.exists(genbank_output_dir): #creates a output path for all mitochondrial nodes
    os.mkdir(genbank_output_dir)
gff_output_dir = "gffs/" # Output for all  blast results
if not os.path.exists(gff_output_dir):
    os.makedirs(gff_output_dir)



for line in open(path_to_file, 'r'):
    if line[0] == "#":
        continue
    else:
        species_name, Group, SubGroup, Type, RefSeq, accession, Size_Kb, GC, Protein, rRNA, tRNA, Other_RNA, Gene, Pseudogene, Release_date, modify_date = line.split("\t")
        species_name = species_name.replace(" ", "_")
        print (species_name, accession)
        command = ['curl', 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + accession + '&rettype=gb&retmode=text', '-o', genbank_output_dir + species_name.replace(".", "") + "_" + str(count) + '.gbk']
        fasta_command = ['curl', 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + accession + '&rettype=fasta&retmode=text', '-o', fasta_output_dir + species_name.replace(".", "") + "_" + str(count) + '.fasta']
        sp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sp.communicate()
        sd = subprocess.Popen(fasta_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sd.communicate()
        count += 1

# def extract_file_names(search_term, path_to_dir):
#     '''extracts a set of file names that contain a search term and returns it as a list'''
#     filenames = []
#     file_list = os.listdir(path_to_dir)
#     for file_name in file_list:
#         if search_term in file_name:
#            filenames.append(file_name)
#     return filenames
#
# gb_files = extract_file_names(".gbk", "genbank/")
#
# for file in gb_files:
#     print (file)
#     genbank_to_gff_command = ['genbank_to_gff.py', file]
#     std = subprocess.Popen(genbank_to_gff_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     std.communicate()
