#!/usr/bin/python3

#author --> Joseph Sevigny and Nate Ennis
#purpose --> run the art command to simulate a metagenomic assembly
        # --> art_illumina -sam -i Sequence_genomic.fna -p -l 250 -f 10 -m 650 -s 150 -o Sequence_l-250_f-10_m-650_s-150
#dependencies 


# input file is a space seperated list of genbank_accessions, wanted_coverage, and sample_name
#example_file
#NC_007795 10 Staphylococcus_aureus

import sys, os, subprocess

art_flag = True

if len(sys.argv) > 2: #check for flag to run art
    if sys.argv[2] == 'no-art':
        art_flag = False


fasta_dir = 'fasta_dir/' #sys.argv[1]
if not os.path.exists(fasta_dir): #creates a output path for all mitochondrial nodes
    os.mkdir(fasta_dir)

if art_flag:
    output_dir = 'art_output/'
    if not os.path.exists(output_dir): #creates a output path for all mitochondrial nodes
        os.mkdir(output_dir)


if art_flag:
    output_log = open('art_log.txt', 'w')

art = "/home/mcbs913kd/art_bin_CC_Cake/art_illumina"
count = 0

input_file = open(sys.argv[1], 'r')

for line in input_file.readlines(): # goes through our file one line at a time
    count += 1
    print ("\n\n#####Genome ", count, "#####\n")
    line = line.rstrip() # removes new line character form the line
    if art_flag:
        genome, coverage, name = line.split(' ')
    else:
        genome = line.split(' ')[0]
        name = 'unknown'
    genome_file = fasta_dir + genome + '.fasta'
    if art_flag:
        output_name = output_dir + genome + '_'


    print ('downloading genome --> ', genome, name)
    fasta_command = ['curl', 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + genome + '&rettype=fasta&retmode=text', '-o', fasta_dir + genome + '.fasta']
    gff_command = ['curl', 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + genome + '&rettype=fasta&retmode=text', '-o', fasta_dir + genome + '.gff']
    sd = subprocess.Popen(fasta_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sd.communicate()
    sd = subprocess.Popen(gff_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sd.communicate()
    
    print ('genome_download_successful:', fasta_dir + genome + '.fasta\n')

    if art_flag:
        art_command = [art, "-sam", "-i", genome_file, "-p", "-l", "250", "-f", coverage, "-m", "650", "-s", "150", "-o", output_name]
        print ("running art command on {}: {}: coverage ({})".format(genome,name, coverage))
        output_log.writelines('genome\t' + genome_file + '\tcoverage\t' + coverage + '\t' + name + '\n' + str(art_command) + '\n')
        sp = subprocess.Popen(art_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sp.communicate() # runs one at a time instead of all at once!!!!!
        #could improve the script by making it multi threaded here

print ("\n\nAll done, have a great day!! enjoy this random thing!!!\n\n")

os.system('random')



