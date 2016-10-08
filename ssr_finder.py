#!/usr/bin/python3

import os, sys, re
from Bio import SeqIO

path_to_misas = "NEM_misas/edit_misas/"
path_to_contigs = "NEM_contigs/"
path_to_sams = "NEM_sams/"
outputfile = open("NEW_SSR_OUTPUT", 'w')

common_thresh = 8

outputfile.writelines("ID\tFOREIGN_ID\tUPSTREAM\tDOWNSTREAM\tSSR\n")

def extract_file_names(search_term, path_to_dir):
    '''extracts a set of file names that contain a search term and returns it as a list'''
    filenames = []
    file_list = os.listdir(path_to_dir)
    for file_name in file_list:
        if search_term in file_name:
           filenames.append(path_to_dir + file_name)
    return filenames

misa_list = extract_file_names(".misa", path_to_misas)
contig_list = extract_file_names(".fasta", path_to_contigs)
sam_list = extract_file_names(".sam", path_to_sams)

print (misa_list, contig_list, sam_list)

contigs_misaNsam = {}
for contig in contig_list:
    tmp, finder_tmp = contig.split("/"); finder, tmp = finder_tmp.split(".")
    for misa in misa_list:
        if finder in misa and finder in contig:
            contigs_misaNsam[contig] = [misa]
# for sam in sam_list:
#     tmp, finder_tmp = sam.split("/"); finder, tmp = finder_tmp.split(".")
#     for contig in contig_list:
#         for misa in misa_list:
#             if finder in misa and finder in contig and finder in sam:
#                 contigs_misaNsam[contig] = [misa,sam]
print ("\n",contigs_misaNsam)

def misa_extractor(line):
    id, ssr_nr, ssr_type, ssr, size, start, end = line.split("\t")
    return (id, ssr_type, ssr, int(start), int(end))

def rev_comp(seq):
    """Reverses, complements and returns sequence"""
    rev_seq = seq[::-1]
    compliment_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    rev_comp_seq = ''
    for nuc in rev_seq:
        if nuc in ['A', 'T', 'G', 'C']:
            rev_comp_seq += compliment_dict[nuc]
    return rev_comp_seq

def extract_region(fasta, header_search, start, stop):
    """ PURPOSE: Extracts a specified region from a fasta using header_search, coordinates, and direction(+-),
                        and writes to append able output
        USAGE: extract_region(fasta, header_search, start, stop, output_name, direction, header_addition)"""
    for seq_record in SeqIO.parse(fasta, "fasta"):
        if header_search in str(seq_record.id):
            my_seq = (str(seq_record.seq[start-1:stop]))
    begin_flank = my_seq[:30]; end_flank = my_seq[-30:]; ssr = my_seq[30:-30]
    return (begin_flank, end_flank, ssr)

# def locate_ssrs(old_contig,contig_file, begin_flank, end_flank, count):
#     for seq_record in SeqIO.parse(contig_file, "fasta"):
#         #str_expression = begin_flank + "(w/{1,500})" + end_flank
#         str_expression = begin_flank + "(\w*)" + end_flank
#         rev_begin = rev_comp(begin_flank); rev_end = rev_comp(end_flank)
#         rev_str_expression = rev_begin + "(\w*)" + rev_end
#         print (rev_str_expression)
#         print (str_expression)
#         if begin_flank in str(seq_record.seq) and end_flank in str(seq_record.seq):
#             reg_exp = re.compile(str_expression)
#             match = reg_exp.search(str(seq_record.seq)) #findall maybe
#             print ("found it!!!! ---> ", match)
#             print (seq_record.id)
#             count += 1
#             new_line = "\n" + old_contig + "\t" + contig_file + seq_record.id + "\t" + begin_flank + "\t" + end_flank + "\t" +  str(match)
#         elif rev_begin in str(seq_record.seq) and rev_end in str(seq_record.seq):
#             rev_reg_exp = re.compile(rev_str_expression)
#             match = rev_reg_exp.search(str(seq_record.seq))
#             print ("found reverse!!!! ---> ", match)
#             count +=1
#             new_line = "\n" + old_contig + "\t" + contig_file + seq_record.id +  "\t" + begin_flank + "\t" + end_flank + "\t" +  str(match)
#         return new_line, count

flank = 30
for current_contig, misaNsam in contigs_misaNsam.items():
    with open(misaNsam[0], 'r') as m:
        for line in m:
            if line.startswith("ID"):
                continue
            else:
                current_block = ""
                id, ssr_type, ssr_id, start, end = misa_extractor(line)
                start = start - flank; end = end+flank
                begin_flank, end_flank, ssr = extract_region(current_contig,id,start,end)
                tmp, old = current_contig.split("/")
                current_block += old + "\t" + old + "\t" + begin_flank + "\t" + end_flank + "\t" + ssr
                count = 0
                if begin_flank:
                    for my_contig in contig_list:
                        #new_line, count = locate_ssrs(current_contig,my_contig,begin_flank,end_flank,count)
                        for seq_record in SeqIO.parse(my_contig, "fasta"):
                            if my_contig == current_contig:
                                continue
                            else:
                                #str_expression = begin_flank + "(w/{1,500})" + end_flank
                                str_expression = begin_flank + "(\w{1,300})" + end_flank
                                rev_begin = rev_comp(begin_flank); rev_end = rev_comp(end_flank)
                                rev_str_expression =  rev_end + "(\w{1,300})" + rev_begin
                                # print (rev_str_expression)
                                # print (str_expression)
                                if begin_flank in str(seq_record.seq) and end_flank in str(seq_record.seq):
                                    reg_exp = re.compile(str_expression)
                                    match = reg_exp.findall(str(seq_record.seq)) #findall maybe
                                    #print ("found it!!!! ---> ", match)
                                    count += 1
                                    match = str(match)
                                    tmp, old = current_contig.split("/")
                                    tmp, new = my_contig.split("/")
                                    new_line = "\n" + old + "\t" + new + "\t" + begin_flank + "\t" + end_flank + "\t" +  match
                                    current_block += new_line
                                elif rev_begin in str(seq_record.seq) and rev_end in str(seq_record.seq):
                                    rev_reg_exp = re.compile(rev_str_expression)
                                    match = rev_reg_exp.findall(str(seq_record.seq))
                                    #print ("found reverse!!!! ---> ", match)
                                    match = rev_comp(str(match))
                                    #print (match)
                                    count +=1
                                    tmp, old = current_contig.split("/")
                                    tmp, new = my_contig.split("/")
                                    new_line = "\n" + old + "\t" + new + "\t" + begin_flank + "\t" + end_flank + "\t" +  match
                                    current_block += new_line
                    if count >= common_thresh:
                        print (current_block)
                        outputfile.writelines(current_block + "\n" + "\n")

def do_blast(command_list):
    '''performs blasts and returns the results'''
    sp = subprocess.Popen(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    blast = sp.communicate()
    blast = str(blast)
    return (blast)
