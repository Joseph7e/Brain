#!/usr/bin/python3
#
#  fastaRead.py fastafilec
#               Read a FastA file representing genome sequence data,
#                   read 1 sequence header and its data at a time
#                   do <something> for the sequence
#                   output <something>
#               This version computes the length of the sequence and outputs
#               the sequence header line with the length info appended.
#   rdb
#   01/07/15
#
#   01/08/15 Upgraded to Python 3.4:
#
import os    # module that handles OS features
import sys   # sys module needed to access command line information

#----------- global variables -----------------------
usageMsg = '''

Usage: fastaRead fastafile
          Extract each sequence from a fastafile into a header string and
          a sequence string.
          <do something to the sequence>
          This version computes its length and prints the seqId and length
          to standard output.

'''

#-------------------------- usage() ---------------------------------------
#----- usage test 
#      If we have no arguments or the first argument is -h
#      print some "usage" information and quit.
#
# Note that the sys.argv array includes the "programName" in positions 0
#      and command line arguments begin at position 1.

def usage():
    if len( sys.argv ) < 2 or sys.argv[ 1 ] == "-h":  # command 
        print( usageMsg )
        exit( 0 )
#---------------------------------------------------------------------    
#----------------------- processSequence -----------------------------
# Given a header and a sequence as a single character string with no line feeds,
#   do whatever processing is desired for that sequence.
# The starting code just computes the sequence size and prints the seq id and
#   size to standard output.
#
def processSequence( header, sequence ):
    # ------------------------------------------
    #  Supplement or replace the lines below with the 
    #  sequence processing and output you want to do.
    #
    basesCount = len( sequence )

    # extract the sequence id from the header
    header = header + " " # add a space at end; guarantees that index finds one   
    spacePos = str.find( header, ' ' ) # find first space in header
    seqIdEnd = spacePos - 1            # index of end of sequence id
    seqId = header[ 1 : seqIdEnd ]     # [ start : end ] is python "slice" operator
                                       # slice works for both strings and arrays.
    #Task 1) Calculate GC-content
    sequence = sequence.upper() # Insures all uppercase bps
    g_count = sequence.count("G"); c_count = sequence.count("C")
    t_count = sequence.count("T"); a_count = sequence.count("A") # counts each nucleotide individually
    length_of_standard_nucs = (g_count + c_count + t_count + a_count)
    GC_content = (g_count + c_count) / length_of_standard_nucs

    #Task 2) calculate number and percentages of non gcta nucleotides
    number_of_non_standard = basesCount - length_of_standard_nucs
    percentage_non_standard = number_of_non_standard/basesCount*100

    # ----- alternative approach
    # nuc_list = ["A","C","T","G"]
    # count = 0
    # for nuc in sequence:
    #     if nuc not in nuc_list:
    #         count += 1
    # percentage_non_standard = count/basesCount*100

    # print automatically adds newline appropriate for current OS
    print( 'seqId=' + seqId  + ' | length=' + str( basesCount ) + ' | GC-content = ' + str(round(GC_content,3)) + ' | number of non-ATCG nucs = ' + str(number_of_non_standard) + '(' + str(round(percentage_non_standard,3)) + '%)' )
    return header, basesCount, GC_content, number_of_non_standard, percentage_non_standard
#-----------------------------------------------------------------------
#--------------------------- main --------------------------------------
# check argument usage correctness
usage()

seqFile = sys.argv[ 1 ]

#----- open file
try:
    inFile = open( seqFile, 'r' )
except ( OSError, IOError ) as e:
    print( 'Unable to open ', seqFile )

#----- read file and process sequences
# first line better be a sequence header
header = inFile.readline()

if not header.startswith( '>' ):
    print( "*** ERROR: Expecting sequence header, got:\n{0}".format( header ),
           end="", file=sys.stderr )
    print( "****** is this a fasta file??? ", file=sys.stderr ) 
    sys.exit()

#----- read and process the sequences until done
number_of_sequences = 0 ## the following lists keep track of all sequence statistics, could also use dict
sequence_length_list = []
gc_content_list = []
non_GCAT_number_list = []
non_GCAT_percentage_list = []
sequence_headers = []
while header != '':
    number_of_sequences += 1 # keeps track of number of headers
    header = header.rstrip( os.linesep )   # delete line separator for any OS   
    #-------------------------------------------------
    # Non-regex code for extracting the sequence Id from the header
   
    #
    seq = ''
    line = inFile.readline()  # returns empty string on eof
    
    while  line != '' and not line.startswith( '>' ):
        line = str.rstrip( line ) # delete trailing white space including lf
        seq += line               # append line sequence data to seq
        line = inFile.readline()  # read next line
        
    sequence_header, sequence_length, gc_content, non_GCAT_number, non_GCAT_percentage = processSequence( header, seq )
    sequence_length_list.append(sequence_length)
    gc_content_list.append(gc_content)
    non_GCAT_number_list.append(non_GCAT_number)
    non_GCAT_percentage_list.append(non_GCAT_percentage)
    sequence_headers.append(sequence_header)
    #------------------------------------------
    header = line    # last line read is either next header or null

#-----File Statistics ----
def calculate_average(list):
    total_length = 0
    for l in list:
        total_length+=l
    average_length = total_length/len(list)
    return average_length

def calculate_standard_deviation(list, average):
    tmp_standard_list = []
    for l in list: #for each number: subtract the Mean and square the result.
        tmp_standard_list.append((l-average)**2)
    tmp_standard_length = 0 #work out the mean of those squared differences.
    for l in tmp_standard_list:
        tmp_standard_length+=l
    tmp_mean_standard = tmp_standard_length/number_of_sequences
    standard_deviation = tmp_mean_standard**(0.5) # Take the square root to get standard deviation
    return standard_deviation

# calculate average sequence length
avg_length = calculate_average(sequence_length_list)
#calculate standard deviation of sequence length
sd_length = calculate_standard_deviation(sequence_length_list,avg_length)

#calculate average gc-content
avg_gc = calculate_average(gc_content_list)
#calculate sd of gc_content
sd_gc = calculate_standard_deviation(gc_content_list,avg_gc)

#Add all non_gcat characters
total_non_gcat = 0
max_non_gcat = 0
max_sequence_header = ''
max_sequence_length = 0
count = 0
total_number_of_sequences_w_non_GCAT = 0
for n in non_GCAT_number_list: # keep track of non-gcat characters
    count += 1
    total_non_gcat+= n
    if n >= 1:
        total_number_of_sequences_w_non_GCAT += 1
    if n > max_non_gcat:
        max_non_gcat = n
        max_sequence_header = count #pair the sequence header with the max gcat number
tmp_count = 0
for h in sequence_headers:
    tmp_count += 1
    if tmp_count == count:
        max_sequence_header = h
tmp_count = 0
for l in sequence_length_list: #pair the sequence length with max gcat to calculate percentage
    tmp_count += 1
    if tmp_count == count:
        max_sequence_length = l

max_percentage_non_CGAT = (max_non_gcat/max_sequence_length)*100
avg_non_gcat = calculate_average(non_GCAT_number_list)
sd_non_gcat = calculate_standard_deviation(non_GCAT_number_list, avg_non_gcat)
percent_w_non_gcat = (total_number_of_sequences_w_non_GCAT/number_of_sequences)*100

print ("\n\t\t\t-----FASTA STATISTICS-----\nNumber of sequences = {0}".format(number_of_sequences))
print ("Average sequence length --> {0} bp (sd = {1})".format(avg_length, round(sd_length,2)))
print ("Average GC-Content --> {0} (sd = {1})".format(round(avg_gc,2), round(sd_gc,2)))
print ("Total number of non-GCAT characters in all sequences --> {0} (avg = {1}, sd = {2})".format(total_non_gcat, round(avg_non_gcat,2), round(sd_non_gcat,2)))
print ("Percentage of sequences with 1 or more non-GCAT characters --> ", str(round(percent_w_non_gcat,2)) + "%")
print ("Max number of non-GCAT nucleotides in a sequence ({}) --> {} ({}% of seq length)".format(max_sequence_header,max_non_gcat, round(max_percentage_non_CGAT,2)))

print( "\t\t\t-----------------------------------" )
#-------------------- end main -------------------------------------
