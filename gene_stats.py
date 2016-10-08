#!/usr/bin/python3

#takes a directory of mulifastas and outputs a exl file with statistics

import sys, os, re, subprocess
from Bio import SeqIO


###INPUT_STUFF####

fasta_dir = "/home/mcbs913kd/mcbs913sh/downloaded_genomes/fastas/" ##fasta dir used for gene_extractor
input_dir = "/home/genome/joseph7e/gene_16S/analysis_gene_extracting/gene_extract_gyrase/"

output_add = sys.argv[1]

variant_analysis = False


output_dir = output_add + '/'
new_files_dir = output_dir + 'filtered_fastas/'
removed_files_dir = output_dir + 'removed_fastas/'

all_new_fastas_output = open(output_add + "_best.fasta", 'a')
####OUTPUT_STUFF

length_cutoff = 1000


print ("input_fasta -->\t" + input_dir + "\noriginal_fasta_dir --> \t" + fasta_dir + "\noutput_addition -->\t" + output_add)
print ("length_cutoff -->\t" + str(length_cutoff))

if not os.path.exists(output_dir):
    os.mkdir(output_dir)
if not os.path.exists(new_files_dir):
    os.mkdir(new_files_dir)
if not os.path.exists(removed_files_dir):
    os.mkdir(removed_files_dir)

gene_and_length_output = open(output_dir + 'gene_and_length_all.tsv', 'w')
gene_and_length_best_output = open(output_dir + 'gene_and_length_best.tsv', 'w')
variation_table = open(output_dir + 'variation_table_best.tsv', 'w')


####

gene_and_length_output.writelines("sample_id\tfasta_file\tgene_id\tlength\tstart\tstop\n")
variation_table.writelines("sample_id\tcopy_number\tlength_low\tlength_high\t%variation_avg\tvariation_max\n") # how to deal with gaps or incomplete genes

gene_dir_dict = {}

fasta_list = os.listdir(fasta_dir)
gene_fasta_list = os.listdir(input_dir)

for fasta in fasta_list:
    greater_than_count = 0
    less_than_count = 0
    for f in gene_fasta_list:
        if fasta[:13] in f:
            start_stop_list = []
            for seq_record in SeqIO.parse(input_dir + "/" + f, "fasta"):
                #4exon_GCF_000010305.1_ASM1030v1_genomic.fna_1882855_1884385_+1531_NC_012489.1
                my_seq = (str(seq_record.seq))
                sample_id = fasta[:13]
                contig_id = seq_record.id[-11:]
                gene_id = fasta[14:-6]
                # print (seq_record.id)
                # print (sample_id)
                # print (contig_id)
                # print (gene_id)
                start_stop = re.findall(r"_([0-9]*_[0-9]*_.)", seq_record.id)
                start, stop, direction = start_stop[0].split('_')
                header = ("\n>" + sample_id + '_' + contig_id +  '_' + gene_id + "_" + direction + '_start_' + start + 'stop_' + stop + 'length_' + str(len(my_seq)) + "\n")
                #print (start, stop, direction)

                if start_stop in start_stop_list:
                    continue # could print to a tmp file of ignored stuff
                else:
                    start_stop_list.append(start_stop)
                    if len(my_seq) >= length_cutoff:
                        greater_than_count += 1
                        greater_than_output =  open(new_files_dir + sample_id + ".fasta",'a')
                        greater_than_output.writelines(header + str(my_seq))
                        all_new_fastas_output.writelines(header + str(my_seq))
                        gene_and_length_best_output.writelines(sample_id+'\t'+fasta+'\t'+gene_id + '\t' + str(len(my_seq))+ "\t" + start + "\t" + stop + "\n")
                    else:
                        less_than_count += 1
                        less_than_output = open(removed_files_dir +  sample_id, 'a')
                        less_than_output.writelines(header + str(my_seq))
                    gene_and_length_output.writelines(sample_id+'\t'+fasta+'\t'+gene_id + '\t' + str(len(my_seq))+ "\t" + start + "\t" + stop + "\n")
    print (fasta, greater_than_count, less_than_count)


if variant_analysis:
    for f in os.listdir(new_files_dir):
        r = open(new_files_dir + f, 'r')
        copy_number = 0
        variation_list = []
        high_variation = 0
        avg_variation = 0
        length_low = 0
        length_high = 0
        for line in r.readlines():
            if ">" in line: # GCF_000005825_NC_013791.2_exon_+_start_97241stop_98796length_1556
                copy_number+= 1
                lister = re.findall(r"length_(.*)\n", line)
                length = lister[0]
                if int(length) <= length_low or length_low == 0:
                    length_low = int(length)
                if int(length) >= length_high:
                    length_high = int(length)
        if copy_number == 1:
            variation_table.writelines(f+"\t"+str(copy_number)+"\t"+str(length_low)+"\t"+str(length_high)+"\t"+"100"+"\t"+"100" + "\n")
        else:
            command = ['blastn', '-query', new_files_dir + f, '-subject', new_files_dir + f, '-outfmt', '6']
            sp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            results = str(sp.communicate()[0].decode('ascii'))
            for line in results.split("\n"):
                try:
                    query, subject, variation, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9 = line.split("\t")
                except ValueError:
                    continue
                if query == subject:
                    continue
                else:
                    variation_list.append(eval(variation))
                    if eval(variation) <= high_variation or high_variation == 0:
                        high_variation = eval(variation)
            print (variation_list, f )
            try:
                avg_variation = sum(variation_list)/len(variation_list)
            except ZeroDivisionError:
                avg_variation = 0
                print ("missing some blast results")
            variation_table.writelines(f+"\t"+str(copy_number)+"\t"+str(length_low)+"\t"+str(length_high)+"\t"+str(avg_variation)+"\t"+str(high_variation) + "\n")
            if len(variation_list) != (copy_number*copy_number-copy_number):
                print ("missing some blast results")







#print("\nProcess complete, a total of {} contigs in the file greater than {}bp\nScatter plots and tables can be found in {}\nGoodbye".format(len(GC_count_list), length_cutoff, output_dir))


                 # need a regular expression to separate these fields
                ### parse the file to get location length info, a final dataset will group just >1000 bp and give length and variation info

    # for fasta in fasta_list:
    #     if fasta[:-4] in gene_list: # does not work correctly
    #         continue
    #     else:
    #         log.writelines(fasta + "\n")

