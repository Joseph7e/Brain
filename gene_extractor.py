#!/usr/bin/python3

from Bio import SeqIO
import os, sys, gzip, mito_tools, re

##ASSUMPTIONS
# fasta must be unzipped. gff can be either
# MATCHING FILES MUST HAVE IDENTICAL NAMES (extensions = .fasta or .fa or .fna and .gff respectively

###CONFIGS### choose to use sys or config yourself, not both

#sys
# FASTA_dir = sys.argv[2]
# GFF_dir = sys.argv[3]
# output_dir = sys.argv[1]



#configs
FASTA_dir = '/home/mcbs913kd/mcbs913sh/downloaded_genomes/fastas/'
GFF_dir = '/home/mcbs913kd/mcbs913sh/downloaded_genomes/gffs/'
output_dir = 'gene_extract_' + sys.argv[1]

if not os.path.exists(output_dir.replace(" ", "_")):
    os.mkdir(output_dir.replace(" ", "_"))
output_log = open("gene_extract_log_" + sys.argv[1][:-1], 'w')
output_stats = open("gene_stats_" + sys.argv[1][:-1], 'w')




##additional configuartions
flanking = 0 #change to grab flanking regions on either side of the gene


#unique_gene_identifiers = ("16S", "XXXXXXXXXX") # if giving only one be sure to give a blank one as well!!!
#positive_check = ("Rrna", "rrna","rRNA") # ensure this is the gene type
#negative_check = ("CDS", "product=16S", "exon", "region", "tRNA", "gene") # not this gene type

#unique_gene_identifiers = ("cox1", "")
#unique_gene_identifiers = ("cox1", "nad3", "nad1", "cox2", "atp8", "cox3", "nad6", "cob", "atp6", "nad5", "nad4l", "nad4", "nad2")
#unique_gene_identifiers = ("", " ") #will take all genes present
#unique_gene_identifiers = ("trnW", "trnS1", "trnD", "trnC", "trnC", "trnM", "trnL1", "trnL2", "trnY", "trnP", "trnS2", "trnT", "trnH", "trnE", "trnG", "trnK", "trnA", "trnF", "trnQ", "trnR", "trnN", "trnI")
#positive_check = ("CDS", "cds") # ensure these are not in line --> leave as [] for nothing
#unique_gene_identifiers = ("RNA Polymerase Sigma-70", "XXXXXXXX")
#positive_check = ("RNA Polymerase Sigma", "70") # ensure these are not in line --> leave as [] for nothing
#negative_check = ("XXXXXX") # make sure this is not in line --> leave as [] for nothing

### HOUSE-KEEPING GENES ###
#unique_gene_identifiers = ["30S ribosomal protein S12"]
#unique_gene_identifiers = ["30S ribosomal protein S15"]
#unique_gene_identifiers = ["GTPase Der"]
#unique_gene_identifiers = ["ATP synthase subunit delta"]
#unique_gene_identifiers = ["CTP synthase"]
#unique_gene_identifiers = ["DNA gyrase subunit B"] # straight gene identifier
# unique_gene_identifiers = ['translation initiation factor IF-2']
#
#
# positive_check = ("CDS", "gene") # ensure these are not in line --> leave as [] for nothing
# negative_check = ("XXXXXXX","XXXXXXXX") # make sure this is not in line --> leave as [] for nothing
#

unique_gene_identifiers = ['CRISPR spacer']

positive_check = ("ncRNA", "gene") # ensure these are not in line --> leave as [] for nothing
negative_check = ("exon","XXXXXXXX") # make sure this is not in line --> leave as [] for nothing

output_log.writelines("#Fasta Dir -->\t" + FASTA_dir + "\n#GFF Dir -->\t" + GFF_dir + "\n#unique gene ids -->\t" + str(unique_gene_identifiers) + "\n#positive check -->\t" + str(positive_check) + "\n#negative check -->\t" + str(negative_check) + "\n")
output_stats.writelines("#GCF_file\tNode_name\tstart\tstop\tlength\n")

print ("#Fasta Dir -->\t" + FASTA_dir + "\n#GFF Dir -->\t" + GFF_dir + "\n#unique gene ids -->\t" + str(unique_gene_identifiers) + "\n#positive check -->\t" + str(positive_check) + "\n#negative check -->\t" + str(negative_check))

#Obtain gffs and fastas using the extract_file_names definition
gffs = mito_tools.extract_file_names("gff", GFF_dir)
fastas = mito_tools.extract_file_names(".fa", FASTA_dir)
if len(fastas) == 0:
    fastas = mito_tools.extract_file_names(".fna", FASTA_dir)


#create a dictionary of fastas and matching annotations:
fa_gff_dict = mito_tools.match_fasta_and_gff(fastas,gffs) ##dictionary containing all fastas:gffs

print ("number of starting fastas = {}\tnumber of starting gffs = {} \ntotal matching = {}".format(len(fastas), len(gffs), len(fa_gff_dict)))

output_log.writelines("#number of starting fastas = {}\tnumber of starting gffs = {} \n#total matching = {}".format(len(fastas), len(gffs), len(fa_gff_dict)))

output_set = set() #what is this for????

for current_fasta, current_gff in fa_gff_dict.items():
    current_sample = current_fasta[:15]
    gene_count = 0 #in case more than one gene is enetered at a time
    # current_gene_file = output_dir + file_name + unique_gene_identifiers[0] + ".fasta"
    output_log.writelines(("#\n#\n#" + current_fasta + "\t" + current_gff))

    if current_gff.endswith(".gz"):
        g = gzip.open(GFF_dir+current_gff, 'rt') #checks for zipped gff file
    else:
        g = open(GFF_dir + current_gff, 'r')

    print (mito_tools.mess_with_font.Green + "\nWorking on details for each gene in the file --> " + mito_tools.mess_with_font.ENDC, current_fasta)
    for line in g.readlines():
        if line[0] == "#":
            continue
        else:
            for gene in unique_gene_identifiers:
                negative_flag = "off"
                positive_flag = "off"
                if gene in line:
                    header, program, gene_identity, start, stop, tmp, direction, tmp2, product = line.split("\t")
                    output_set.add(gene_identity)
                    #you would need to keep track of each of the unique gene_ids and add it to corrosponding list
                    current_gene_file = output_dir + "/" +  current_fasta[:13] + '_' + gene.replace(" ", "_") + '_' + gene_identity + ".fasta"
                    if negative_check:
                        for n in negative_check:
                            if n in gene_identity:
                                negative_flag = "on"
                            else:
                                continue
                    else:
                        negative_flag = "off"
                        continue

                    if positive_check:
                        for p in positive_check:
                            if p in gene_identity:
                                positive_flag = "on"
                            else:
                                continue
                    else:
                        positive_flag = "on"
                        continue

                    if positive_flag == "on" and negative_flag == "off":
                        gene_count += 1
                        output_log.writelines(line)
                        start = int(start); stop = int(stop)
                        output_stats.writelines(current_sample + '\t' + header + '\t' + gene_identity + '\t' + str(start) + '\t' + str(stop) + '\t' + direction + '\t' + str(stop - start) + "\n")
                        mito_tools.extract_region(FASTA_dir + current_fasta, header, start-flanking, stop+flanking, current_gene_file, direction, str(gene_count)+ "_" + gene_identity+ "_" + current_fasta+"_"+str(start)+"_"+str(stop)+"_" +direction)
                        print (gene, header, gene_identity, str(start), str(stop), direction)
                    else:
                        output_log.writelines("filtered_out__" + line)
                        #print (header, gene_identity, start, stop, direction, mito_tools.mess_with_font.Red + str(negative_check) + mito_tools.mess_with_font.ENDC)


print (output_set)
