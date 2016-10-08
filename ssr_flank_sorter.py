#!/usr/bin/python3
import subprocess, sys, os, re

#INPUT from flank finder and line numbers added to each line with a tab

input_file = "/home/genome/joseph7e/SSRs/new_analyses_march/nemertea_microsats_with_numbers_no_spaces"

file_with_revs = input_file + "_and_revs.txt"

output_file = open(file_with_revs, 'w')

os.mkdir("all_matches/")

agrep_spec = '4'
output_add = "nemertea"


#### ADD A LIST OF HEADER IDENTIFIERS HERE AND TURN INTO A DICTIONARY BELOW

header_tag = "NEM-"
header_ids = {'A02_N50':0, 'A12_N50':0, 'A3_N50':0, 'A4_N50':0, 'A6_N50':0, 'A7_N50':0, 'B02_N50':0, 'B03_N50':0, 'B1_N50':0, 'B7_N50':0, 'E9_N50':0, 'G2_N50':0, 'H2_N50':0}
atlantic = ['A02_N50', 'A12_N50', 'A3_N50', 'A6_N50', 'A7_N50', 'B02_N50', 'B03_N50', 'B1_N50']
pacific = ['A4_N50', 'B7_N50', 'E9_N50', 'G2_N50', 'H2_N50']

all_labels = atlantic+pacific


#### MAKE OUTPUTS TELL US SOMETHING ABOUT EACH GROUPING, ADD CONDITIONS ABOUT MATCHES TO make a bunch of outputs

OUTPUT1 = open(output_add + "_greater_than_half.fs", 'w')
OUTPUT2 = open(output_add + "_all_samples.fs", 'w')
OUTPUT3 = open(output_add + "_most_samples.fs", 'w')
OUTPUT4 = open(output_add + "_greater_than_half_p_and_a.fs", 'w')
OUTPUT5 = open(output_add + "_high_copy_number_matches.fs", 'w')
OUTPUT6 = open(output_add + "_atlantic.fs", 'w')
OUTPUT7 = open(output_add + "_pacific.fs", 'w')


all_matches_output = open(output_add + "_all_matches.fs", 'w')

# output = open("perfect_matches4.txt", 'w')
# output2 = open("notsoperfect_matches4.txt", 'w')
# output3 = open("most_samples_matches4.txt", 'w')

def rev_comp(seq):
    """Reverses, complements and returns sequence"""
    rev_seq = seq[::-1]
    compliment_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    rev_comp_seq = ''
    for nuc in rev_seq:
        if nuc in ['A', 'T', 'G', 'C']:
            rev_comp_seq += compliment_dict[nuc]
    return rev_comp_seq



with open(input_file,) as i:
    for line in i:
        try:
            z,a,b,c,d,e,f,g,front,back,ssr = line.split("\t")
            output_file.writelines(line.strip() + "\t" + rev_comp(back) + "\t" + rev_comp(front) + "\t" + rev_comp(ssr) + "\n")
        except ValueError:
            #print (line)
            continue

output_file.close()

growing_set = set()

with open(file_with_revs) as file:
    print ("Determining reverse complements for flanking regions and ssr")
    reg = ('\d*')
    for line in file:
        counter = re.match(reg, line).group()
        if counter in growing_set:# or line.startswith("SAMPLE"):
            print ("line skipped --> It has been found before")
            continue
        #print(line)
        else:
            ###HERE ADD THE SAMPLE NAMES AND FIX BELOW######
            ####
            #A02_N50 = 0; A12_N50 = 0; A3_N50 = 0; A4_N50 = 0; A6_N50 = 0; A7_N50 = 0; B02_N50 = 0; B03_N50 = 0; B1_N50 = 0; B7_N50 = 0; E9_N50 = 0; G2_N50 = 0; H2_N50 = 0

            #B4_N50 = 0; D8_N50 = 0; G1_N50 = 0; C11_N50 = 0; E7_N50 = 0; C12_N50 = 0; E8_N50 = 0; C2_N50 = 0; E9_N50 = 0
            header_ids = {'A02_N50':0, 'A12_N50':0, 'A3_N50':0, 'A4_N50':0, 'A6_N50':0, 'A7_N50':0, 'B02_N50':0, 'B03_N50':0, 'B1_N50':0, 'B7_N50':0, 'E9_N50':0, 'G2_N50':0, 'H2_N50':0}
            file_name = ''
            try:
                count = 0
                atlantic_count = 0
                pacific_count = 0
                z,a,b,c,d,e,f,g,front, back, ssr, rev_front, rev_back, rev_ssr = line.split("\t")
                grepper = front+"\t"+back
                command = ['agrep', '-'+agrep_spec, grepper, file_with_revs]

                sp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                grep_results = str(sp.communicate()[0].decode('ascii'))
                #count_NODES = grep_results.count("NODE")
                for k, v in header_ids.items():
                    g_count = grep_results.count(header_tag + k)
                    if g_count:
                        header_ids[k] += 1
                        count += 1
                        if k in pacific:
                            pacific_count += 1
                        elif k in atlantic:
                            atlantic_count += 1
                count_NODES = 0
                for text in grep_results.split("\n"):
                    count_NODES+=1
                    new_count = re.match(reg, text).group()
                    growing_set.add(new_count)
                count_NODES -= 1
                print (counter, "Total matching nodes --> ", count_NODES, " . . . Unique Samples --> ", count)
                if count >= ((len(header_ids)/2)-2) and count_NODES <= count+2:
                    print (grep_results)
                    OUTPUT1.writelines(grep_results + "#\n#\n")# now just change output files
                    if pacific_count > 0 and atlantic_count > 0:
                        OUTPUT4.writelines(grep_results + "#\n")
                elif count >= 2 and count_NODES > count+2:
                    OUTPUT5.writelines(grep_results + "#\n")
                if count >= 8 and count_NODES <= count+2:
                    OUTPUT3.writelines(grep_results + "#\n")
                if count >= len(header_ids) and count_NODES <= count+2:
                    OUTPUT2.writelines(grep_results + "#\n")
                if atlantic_count >= (len(atlantic)-2) and pacific_count == 0 and count_NODES <= count+2:
                    OUTPUT6.writelines(grep_results + "#\n")
                if pacific_count >= (len(pacific)-4) and atlantic_count == 0 and count_NODES <= count+2:
                    OUTPUT7.writelines(grep_results + "#\n")
                if count >= 2 and count_NODES <= count+1:
                    all_matches_output.writelines(grep_results + "#\n")

                    for s in all_labels:
                        value = header_ids[s]
                        if value > 0:
                            n = s.replace("_N50", '')
                            file_name += n + "_"
                    with open("all_matches/" + file_name + 'matches', 'a') as m:
                        m.writelines(grep_results + "#\n")
                    print (file_name + 'matches')



            except ValueError:
                continue

# OUTPUT1 = open(output_add + "_greater_than_half.fs", 'w')
# OUTPUT2 = open(output_add + "_all_samples.fs", 'w')
# OUTPUT3 = open(output_add + "_most_samples.fs", 'w')
# OUTPUT4 = open(output_add + "_greater_than_half_p_and_a.fs", 'w')
# OUTPUT5 = open(output_add + "_high_copy_number_matches.fs", 'w')
# OUTPUT6 = open(output_add + "_atlantic.fs", 'w')
# OUTPUT7 = open(output_add + "_pacific.fs", 'w')





#with open(file_with_revs):

                # B4_N50 = grep_results.count("ANN-B4_N50")
                # if B4_N50:
                #     count += 1
                # D8_N50 = grep_results.count("ANN-D8_N50")
                # if D8_N50:
                #     count += 1
                # G1_N50 = grep_results.count("ANN-G1_N50")
                # if G1_N50:
                #     count += 1
                # C11_N50 = grep_results.count("ANN-C11_N50")
                # if C11_N50:
                #     count += 1
                # E7_N50 = grep_results.count("ANN-E7_N50")
                # if E7_N50:
                #     count += 1
                # C12_N50 = grep_results.count("ANN-C12_N50")
                # if C12_N50:
                #     count += 1
                # E8_N50 = grep_results.count("ANN-E8_N50")
                # if E8_N50:
                #     count += 1
                # C2_N50 = grep_results.count("ANN-C2_N50")
                # if C2_N50:
                #     count += 1
                # E9_N50 = grep_results.count("ANN-E9_N50")
                # if E9_N50:
                #     count += 1
                # A02_N50 = grep_results.count("NEM-A02_N50")
                # if A02_N50:
                #     count += 1
                # A12_N50 = grep_results.count("NEM-A12_N50")
                # if A12_N50:
                #     count += 1
                # A3_N50 = grep_results.count("NEM-A3_N50")
                # if A3_N50:
                #     count += 1
                # A4_N50 = grep_results.count("NEM-A4_N50")
                # if A4_N50:
                #     count += 1
                # A6_N50 = grep_results.count("NEM-A6_N50")
                # if A6_N50:
                #     count +=1
                # A7_N50 = grep_results.count("NEM-A7_N50")
                # if A7_N50:
                #     count +=1
                # B02_N50 = grep_results.count("NEM-B02_N50")
                # if B02_N50:
                #     count +=1
                # B03_N50 = grep_results.count("NEM-B03_N50")
                # if B03_N50:
                #     count +=1
                # B1_N50 = grep_results.count("NEM-B1_N50")
                # if B1_N50:
                #     count +=1
                # B7_N50 = grep_results.count("NEM-B7_N50")
                # if B7_N50:
                #     count +=1
                # E9_N50 = grep_results.count("NEM-E9_N50")
                # if E9_N50:
                #     count +=1
                # G2_N50 = grep_results.count("NEM-G2_N50")
                # if G2_N50:
                #     count +=1
                # H2_N50 = grep_results.count("NEM-H2_N50")
                # if H2_N50:
                #     count +=1


                # if count >= 5 and count_NODES <= count+2:
                #     print (grep_results)
                #     output.writelines(grep_results + "#\n#\n") # now just change output files
                # elif count >= 5 and count_NODES > count:
                #     output2.writelines(grep_results + "#\n")
                # if count >= 8 and count_NODES <= count+2:
                #     output3.writelines(grep_results + "#\n")