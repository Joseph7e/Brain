#!/usr/bin/python3

#Author: Joseph Sevigny
#Affiliation: HCGS
#Purpose: option to run busco on list of samples. max four at a time determine sequence overlap, extract sequences, and align for specified group.
#USAGE: ./busco_dealer.py <sample_commonality> <Path to Contig Directory> <number in common(int)>

import sys, subprocess, os, re
from Bio import SeqIO
#from collections import defaultdict

# sys args will be done with sysargs module
# or excepts any length list seperated by commas. one arg for each side
# dtermines maximum alignments that we can do with atleast a two samples for each,
# done over the whole genome will be significant, greater degree of freedom.
# :))))))))
#
# blast buscos against ncbi???? determine closest reative by grabbing taxonomy
# There must be codes that blobtools used to get taxon id. good test of relatives.
# Determine what the genes do.
#
# Do this with best genes from maker and augustus!! determine closest relatives.
# Are The genes from who we think.
#
# #or make an input config_file
# path/to/contigs \t path/to/busco_run_dir \t group_id

output_dir = "busco_dealer_V3/"
nuc_output_dir = output_dir +  "nucleotides/"
amino_output_dir = output_dir + "amino_acids/"
alignments_dir = output_dir + "amino_alignments/"
nuc_alignment_dir = output_dir + "nucl_alignments/"
#
#
#
if not os.path.exists(output_dir): #creates a directory hierarchy for files
    os.mkdir(output_dir)
if not os.path.exists(nuc_output_dir):
    os.mkdir(nuc_output_dir)
if not os.path.exists(amino_output_dir):
    os.mkdir(amino_output_dir)
if not os.path.exists(alignments_dir):
    os.mkdir(alignments_dir)
if not os.path.exists(nuc_alignment_dir):
    os.mkdir(nuc_alignment_dir)




input_data = open(sys.argv[1], 'r')
sample_dict1 = {}; sample_dict2 = {}; group_paths_dictionary = {} # {'group_id':{'sample':'contig_path','busco_path'}
group_start = ''; group_ids = []

group1 = []
group2 = []

for line in input_data:
    group_id,sample_name,path_to_contigs,path_to_busco = line.rstrip().split(',')
    if group_id not in group_ids:
        group_ids.append(group_id)
    if group_start == '':
        group_start = group_id
    if group_id == group_start:
        group1.append(sample_name)
        sample_dict1[sample_name] = [path_to_contigs,path_to_busco] #{sample_name:[path_to_to_contigs,path_to_busco_runs]}
        group_paths_dictionary[group_id] = sample_dict1
    else:
        group2.append(sample_name)
        sample_dict2[sample_name] = [path_to_contigs,path_to_busco]
        group_paths_dictionary[group_id] = sample_dict2

print ("group1 = {} , samples = {}, totaling {}".format(group_ids[0], group1, len(group1)))
print ("group1 = {} , samples = {}, totaling {}".format(group_ids[1], group2, len(group2)))

#### Amino Acid Variables ####
single_copy_dir = ''




busco_single_copy_amino_dict = {} #busco:[sample_path_to,sample_path_to]

for group,sample_dict in group_paths_dictionary.items():
    for sample, path_list in sample_dict.items():
        contig_path = path_list[0]; busco_path = path_list[1]
        single_copy_dir = busco_path + 'single_copy/'
        current_scg_list = os.listdir(single_copy_dir)
        for scg in current_scg_list:
            if scg in busco_single_copy_amino_dict:
                busco_single_copy_amino_dict[scg].append(group+"|"+sample+"|"+single_copy_dir + scg+"|"+contig_path)
            else:
                busco_single_copy_amino_dict[scg] = []
                busco_single_copy_amino_dict[scg].append(group+"|"+sample+"|"+single_copy_dir + scg+"|"+contig_path)

threshold_count = 0

def muscle_alignment(input_multi_fastas, output_name):
    '''
    example_command = ["muscle", "-clw", "-in", "amino_fasta", "-out", "output_amino"]
    '''
    command = ["muscle", "-clw", "-in", input_multi_fastas, "-out", output_name]
    sp = subprocess.Popen(command)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    alignment = sp.communicate()
def extract_file_names(search_term, path_to_dir):
    '''extracts a set of file names that contain a search term and returns it as a list'''
    busco_run_filenames = []
    file_list = os.listdir(path_to_dir)
    for file_name in file_list:
        if search_term in file_name:
            busco_run_filenames.append(file_name)
    return busco_run_filenames

def rev_comp(seq):
    """Reverses, complements and returns sequence"""
    rev_seq = seq[::-1]
    compliment_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    rev_comp_seq = ''
    for nuc in rev_seq:
        if nuc in ['A', 'T', 'G', 'C']:
            rev_comp_seq += compliment_dict[nuc]
    return rev_comp_seq


for scg, info in busco_single_copy_amino_dict.items():
    #print (len(info))
    print ("#####" + scg + "########" )
    current_list = []
    contig_path_list = []; sample_list = []
    group1_count = 0
    group2_count = 0
    for c in info:
        group,sample,path,contig_path = c.split('|')
        current_list.append(path)
        contig_path_list.append(contig_path)
        sample_list.append(sample)
        if group == group_start:
            group1_count += 1
        else:
            group2_count += 1
    print (group_start, group1_count)
    print ("other_group", group2_count)
    if group1_count >= 2 and group2_count >= 2:
        threshold_count += 1
        output_file = amino_output_dir + scg.replace('.faa', '') + '_'+str(group1_count)+'_'+str(group2_count)+'.fasta'
        output_handle = open(output_file, 'w')
        nuc_output_file = nuc_output_dir + scg.replace('.faa', '') + '_'+str(group1_count)+'_'+str(group2_count)+'.fasta'
        nuc_output_handle = open(nuc_output_file, 'w')
        c_count = 0
        for busco_path in current_list:
            print ('hahahaha busco path', busco_path)
            c_path = contig_path_list[c_count]
            c_sample = sample_list[c_count]
            c_count += 1
            for l in open(busco_path,'r'):
                if l == '\n':
                    continue
                if l.startswith('>'):
                    output_handle.writelines(l)

                    m = l.rstrip()
                    node = re.findall(r'(NODE_[0-9]*_)',m)[0]
                    # start,stop = l.split(':')[3].split('-')
                    print('#',c_sample,c_path)
                    for seq_record in SeqIO.parse(c_path, "fasta"):
                        if node in str(seq_record.id):
                            aug_files = extract_file_names(scg.replace('.fas',''),busco_path[:-31]+'augustus/')
                            direction = ''
                            my_seq = ''
                            first_one = 0
                            last_one = 0
                            for aug in aug_files:
                                seq = ''
                                aug_path = busco_path[:-31]+'augustus/'+aug#.replace(scg,'augustus/'+aug)
                                first_one = 0
                                last_one = 0
                                for line in open(aug_path,'r').readlines():
                                    if line.startswith("#"):
                                        continue
                                    else:

                                        tmp,tmp,type,srt,stp,tmp,d,tmp,tmp = line.split('\t')
                                        if type == "exon":
                                            #print (line)
                                            if first_one == 0:
                                                first_one = srt
                                            if int(stp) > last_one:
                                                last_one = int(stp)
                                            print (aug, "CDS!!!!",srt,stp)
                                            cur_seq = (str(seq_record.seq[int(srt)-1:int(stp)]))
                                            direction = d
                                            seq+=cur_seq
                                if len(seq) > len(my_seq):
                                    my_seq = seq
                                    print ('takeover', aug)
                            print (c_sample, ' length =', len(my_seq),'\n')
                            nuc_output_handle.write('>'+scg+':'+c_sample+':'+str(seq_record.id)+':'+str(first_one)+'-'+str(last_one)+':'+direction+'\n')
                            if direction == "-":
                                rev_seq = rev_comp(my_seq)
                                nuc_output_handle.write(rev_seq)
                            elif direction == "+":
                                nuc_output_handle.write(my_seq)
                            else:
                                print ("error, no direction")
                                sys.exit()
                            nuc_output_handle.write("\n")
                            break
                    #         nuc_output_handle.write(l)
                    #         my_seq = (str(seq_record.seq[int(start)-1:int(stop)]))
                    #
                    #         nuc_output_handle.write("\n")
                    #         # if direction == "-":
                    #         #     rev_seq = rev_comp(my_seq)
                    #         #     my_fasta.write(rev_seq + "\n")
                    #         # else:
                    #         #     my_fasta.write(my_seq + "\n")
                    #         break


                else:
                    output_handle.writelines(l)
        output_handle.close()
        nuc_output_handle.close()
        output_amino = output_file.replace(amino_output_dir,alignments_dir).replace('.fasta','.clw')
        output_nuc = nuc_output_file.replace(nuc_output_dir,nuc_alignment_dir).replace('.fasta','_nuc.clw')
        muscle_alignment(output_file, output_amino)
        muscle_alignment(nuc_output_file, output_nuc)
print (threshold_count)
