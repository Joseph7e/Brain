#!/usr/bin/python3

#AUTHOR: Joseph Sevigny (jlsevigny1@wildcats.unh.edu)
#AFFILIATION: Hubbard Center for Genome Studies _UNH_ _GENETICS_
#PURPOSE: Extract mitochondrial contigs, annotate, and align genes to check annotations
#USAGE: ./mito_contig_finderV2.py <contig(s)> <reference(s)> <output(s)> [option(s)]


import subprocess, sys, os, re
from Bio import SeqIO


help = """
                                #### Welcome to MitoContigFinder ####
        AUTHOR: Joseph Sevigny (jlsevigny1@wildcats.unh.edu) Hubbard Center for Genome Studies
        PURPOSE: Extract mitochondrial contigs, annotate, and align genes to check annotations

        USAGE: ./mito_contig_finderV2.py <contig(s)> <reference(s)> <output(s)> [option(s)]

        contig(s) --> required. path to contig file(s). single argument or comma separated list (assume SPADES fasta)
        reference(s) --> required. path to reference file(s). single argument or comma separated list
        output(s) --> required. name for output fastas. All contigs must have corresponding output
        option(s) --> optional. -h -- display this help menu and exit
                                -v -- verbose mode (default False)
                                -r -- Use the same reference for all files (default False)
                                -a -- runs blast with amino acids for difficult to find mtDNA (default run if it can't find)
                                -M -- run MITOS on found contigs (default False)
                                -G -- extract protein coding genes and create align for each and all (assumes MITOS, default (False)
                                -f -- grab +21bp flanking regions (default = 0) assumes -G
                                -ma -- run muscle alignments for all genes
                                -skip -- skip locating the mitochondrial contigs (default False)
                                -predb

        example usage -->
        ./mito_contig_finderV2.py  files/ANN-B4_contigs.fasta,files/ANN-C2_contigs.fasta files/annelids.fasta ANN-B4,ANN-C2 -v,-G,-f,-r,-ma

        """



#CONFIGS
file_output_dir = "mito_files/"
if not os.path.exists(file_output_dir): #creates a output path for all mitochondrial nodes
    os.mkdir(file_output_dir)
output_dir = "mitochondrial_genomes/"
if not os.path.exists(output_dir): #creates a output path for all mitochondrial nodes
    os.mkdir(output_dir)
blast_output_dir = file_output_dir + "blast_results/" # Output for all  blast results
if not os.path.exists(blast_output_dir):
    os.makedirs(blast_output_dir)
Verbose = False
skipDB = False
Blast = True
Amino = False
MITOS = False
Find_Genes = False
Alignment = False
protein_fasta = "Sc_reference.fasta"
cov_cutoff = 12 #Found to be the best indictor of true mitochondrial node for tblastn results, increase for highly covered nodes

#GENE FINDER STUFF
FASTA_dir = "fastas/"
GFF_dir = "gffs/"




contigs = []; reference = []; output = []; options = [] #Balnk lists of arguments ot be filled with sys.argv

try:
    for option in sys.argv[4].split(","): #optional options
        options.append(option)
except IndexError:
    pass
try:
    for contig in sys.argv[1].split(","): #User comma separated list of contigs or single argument
        contigs.append(contig)
except IndexError:
    print (help) # print help since there are no arguments and exits
    sys.exit()

if "-h" in options or sys.argv[1] == "-h": #print help if needed and exits
    print (help)
    sys.exit()


for ref in sys.argv[2].split(","): #references for blast
    reference.append(ref)
for out in sys.argv[3].split(","): #Must provide output for each contig entered
    output.append(out)

total_files = len(contigs) #keeps track of total input files in case the user enters multiple files
if "-v" in options:
    Verbose = True
if "-r" in options:
    flag = total_files-1
    while flag:
        flag -= 1
        reference.append(reference[0])
if "-a" in options:
    Amino = True
if "-G" in options:
    Find_Genes = True
if "-M" in options:
    MITOS = True
if "-ma" in options:
    Alignment = True
if "-skip" in options:
    Blast = False
if "-preDB" in options:
    skipDB = True


#Create contig dictionary with references and output as keys
contig_dict = {}
count = 0
for contig in contigs:
    contig_dict[contig] = [reference[count], output[count]]
    count += 1


#CLASSES and DEFINITIONS USED THROUGHOUT PROGRAM
class mess_with_font:
            ENDC = '\033[0m'
            BOLD = '\033[1m'
            UNDERLINE = '\033[4m'
            Red = '\033[91m'
            Green = '\033[92m'
            Blue = '\033[94m'
            Yellow = '\033[93m'
            Grey = '\033[90m'
            Default = '\033[99m'

def do_blast(command_list):
    '''performs blasts and returns the results'''
    sp = subprocess.Popen(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #blast = sp.communicate()
    #blast = str(blast)
    blast = str(sp.communicate()[0].decode('ascii'))
    print (command_list)
    print (blast)
    return (blast)

def extract_file_names(search_term, path_to_dir):
    '''extracts a set of file names that contain a search term and returns it as a list'''
    filenames = []
    file_list = os.listdir(path_to_dir)
    for file_name in file_list:
        if search_term in file_name:
           filenames.append(file_name)
    return filenames


def rev_comp(seq):
    """Reverses, complements and returns sequence"""
    rev_seq = seq[::-1]
    compliment_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    rev_comp_seq = ''
    for nuc in rev_seq:
        if nuc in ['A', 'T', 'G', 'C']:
            rev_comp_seq += compliment_dict[nuc]
    return rev_comp_seq

def muscle_alignment(input_multi_fastas, output_name):
    '''
    example_command = ["muscle", "-clw", "-in", "amino_fasta", "-out", "output_amino"]
    '''
    command = ["muscle", "-clw", "-in", input_multi_fastas, "-out", output_name]
    sp = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sp.communicate()

##### RUN THE PROGRAM!!!!
if Verbose:
    print("\n\n\t\t\t{}******* Welcome to Mitochondrial Contig Finder *******\n\n This tool is used to extract mitochondrial genomes from contiguous sequences and align protein coding genes{}.\n".format(mess_with_font.BOLD + mess_with_font.Blue, mess_with_font.ENDC),"=" * 108)
    print("Verbose = ", Verbose, "\nFind Mito Contigs = ", Blast, "\nRun MITOS = ", MITOS, "\nRun tblastn = ",Amino,"\nExtract Protein Coding Genes = ", Find_Genes, "\nMuscle Alignment = ", Alignment)
    print("\nYou have entered a total of {} contigs\n".format(mess_with_font.BOLD + mess_with_font.Red + str(len(contigs)) + mess_with_font.ENDC))

if Blast:
    for contig, details  in contig_dict.items():
        ####Determine which contig(s) contain mitochondrial genomes####
        reference = details[0]
        output_name = details[1]
        blastn_results_file = blast_output_dir + output_name + "_blastn_results" #output for blastn and mito files
        mito_contigs_results_file = output_dir + output_name+ "_mito_contigs"
        writeable_mito_contig_file = open(mito_contigs_results_file, 'w')

        db_command = ["makeblastdb", "-in", contig, "-dbtype", "nucl", "-out", contig + "_db"] #command to construct database
        blast_command = ["blastn", "-query", reference, "-db", contig+"_db", "-outfmt", "6"] #command to complete blastn
        print("{}Current Parameters{}:\nContig File = {}\nReference = {}\nOutput Name = {}".format(mess_with_font.UNDERLINE + mess_with_font.Green, mess_with_font.ENDC, mess_with_font.Red + contig + mess_with_font.ENDC,mess_with_font.Red + reference + mess_with_font.ENDC, mess_with_font.Red + output_name + mess_with_font.ENDC))

        if skipDB:
            continue
        else:
            if Verbose:
                print("\n\nCreating blast database from your contig file . . . DB files will be stored in contig file directory")
            make_db = subprocess.Popen(db_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE) #run database command
            make_db.communicate() #wait till process is complete

        writeable_blast_file = open(blastn_results_file, 'w')
        if Verbose:
            print("\nRunning blastn command . . .", mess_with_font.Red + str(blast_command) + mess_with_font.ENDC)
        blast_results = do_blast(blast_command)
        writeable_blast_file.writelines(blast_results)

        ### DETERMINE UNIQUE MITOCHONDRIAL CONTIGS FROM BLAST RESULTS ####
        my_nodes = []    #locates nodes from blast file
        my_search = re.findall(r"(NODE_\w*.\w*)", blast_results)
        for node in my_search:
            if node not in my_nodes:
                my_nodes.append(node)

        #### DETERMINE COMBINED LENGTH OF ALL NODES ####
        node_length = 0
        my_nodes_list = '' #turns list into a string for searching
        for node in my_nodes:
            my_nodes_list += node
        node_length_list = re.findall(r"th_(\w*)_c", my_nodes_list) #locate the length of the contig from header
        for i in node_length_list: #add up the length of all nodes
            i = int(i)
            node_length += i

        if Verbose:
            print("\nYour mitochondrial nodes are . . . ", mess_with_font.Blue + mess_with_font.BOLD + str(my_nodes) + mess_with_font.ENDC)
            print("The combined length of your nodes(s) = ", mess_with_font.Blue + mess_with_font.BOLD + str(node_length) + mess_with_font.ENDC)

        #### RUN ADDITIONAL tBLASTn IF WANTED OR NECESSARY ####
        if Amino or node_length < 13000:
            tblastn_results_file = blast_output_dir + output_name + "_tblastn_results"
            tblast_command = ["tblastn", "-db", contig +"_db", "-query", protein_fasta, "-db_gencode", "5", "-outfmt", "6", "-out", tblastn_results_file]

            if Verbose:
                if node_length < 13000:
                    print("\nYour node length is smaller than expected")
                print("\n" + mess_with_font.UNDERLINE + mess_with_font.Green + "Running tblastn to pick up additional nodes that may have been missed" + mess_with_font.ENDC)
                print ("Typically mitochondrial will have higher coverage, lets take advantage: Cutoff_coverage =", mess_with_font.Blue + str(cov_cutoff) + mess_with_font.ENDC, "\n")
            if Verbose:
                print("{}Current Parameters{}:\nContig File = {}\nReference = {}\nOutput Name = {}".format(mess_with_font.UNDERLINE + mess_with_font.Green, mess_with_font.ENDC, mess_with_font.Red + contig + mess_with_font.ENDC,mess_with_font.Red + protein_fasta + mess_with_font.ENDC, mess_with_font.Red + output_name + mess_with_font.ENDC))
                print("\nRunning tblastn command . . .", mess_with_font.Red + str(tblast_command) + mess_with_font.ENDC)

            #### SAME AS ABOVE, LOCATING UNIQUE MITO CONTIGS AND ADDING UP LENGTH ####
            tblast_results = do_blast(tblast_command)
            tblast_resultsx = open(tblastn_results_file, 'r')
            for line in tblast_resultsx.readlines():
                cov = re.findall(r"cov_(\w*.\w*)_ID", line)
                for current in cov:
                    current = float(current)
                    if current >= cov_cutoff:
                        node_name = re.findall(r"(NODE_\w*.\w*)", line)
                        for current_node in node_name:
                            if node_name not in my_nodes:
                                my_nodes.append(current_node)

            my_nodes_string = '' #turns list into string for searching
            for node in my_nodes:
                if node not in my_nodes_string:
                    my_nodes_string += node
            node_length_list_again = re.findall(r"th_(\w*)_c", my_nodes_string)
            node_length = 0
            for i in node_length_list_again:
                i = int(i)
                node_length += i

            if Verbose:
                print("...\n...\nYour mitochondrial nodes are . . . ", mess_with_font.Blue + mess_with_font.BOLD + str(my_nodes) + mess_with_font.ENDC)
                print("The combined length of your nodes = ", mess_with_font.Blue + mess_with_font.BOLD + str(node_length) + mess_with_font.ENDC)
                print("...\nExtracting mitochondrial nodes from contig file and writing to output")


        #### EXTRACT MITOCHONDRIAL CONTIG(S) FROM FILE AND RIGHT TO OUTPUT ####
        writeable_mito_contig_file = mito_contigs_results_file
        for current_node in my_nodes:
            with open (writeable_mito_contig_file, 'w') as cfo:
                for seq_record in SeqIO.parse(contig, "fasta"):
                    if current_node in str(seq_record.id):
                        cfo.write(">" + contig + str(seq_record.id) + "\n")
                        cfo.write(str(seq_record.seq) + "\n")

        total_files -= 1
        if Verbose:
            if total_files == 0:
                print ("\n\nProcess complete\n")
            else:
                print ("\n\nCurrent process complete\nMoving to next file\nFiles remaining . . . ".format(total_files), "\n")

#### RUN MITOS ANNOTTAION PROGRAM ####

if MITOS:
    mitos_output_dir = "mitochondrial_genomes/gffs/" #directory to store mitos output
    if not os.path.exists(mitos_output_dir):
        os.mkdir(mitos_output_dir)
    if Verbose:
        print ("Running Mitos to annotate mitochondrial fastas\nThis will take at least an hour")
    fastas = extract_file_names(".fasta", output_dir) #extract the filenames of the mitochondrial genomes
    Ready = False
    if Ready: #Communicating with creator of mitos to get the program to run, he is aware of the error
        for fasta in fastas:
            mitos_command = ["python", "-i", fasta, "-c", "5", "-o", mitos_output_dir]
            sp = subprocess.Popen(mitos_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            sp.communicate()
    else:
        print ("\n {} MITOS COMING SOOON!!!!! {} For now we will work with pre-prepared gffs".format(mess_with_font.Red, mess_with_font.ENDC))


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
            else:
                my_fasta.write(">" + header_addition + str(seq_record.id) + "\n")
                my_seq = (str(seq_record.seq[start-1:stop]))
                if direction == "-":
                    rev_seq = rev_comp(my_seq)
                    my_fasta.write(rev_seq + "\n")
                else:
                    my_fasta.write(my_seq + "\n")

def parse_gff(gff_line):
    '''input a line from a gff and return the important attributes'''
    try:
        header, program, gene_identity, start, stop, tmp, direction, tmp2 = gff_line.split("\t")
    except ValueError:
        header, program, gene_id, start, stop, tmp, direction, tmp2, gene_identitytmp = gff_line.split("\t")
        gene_identity = gene_identitytmp[5:11]
    gene_identity = gene_identity.lower(); gene_identity = gene_identity.replace("nd", "nad").replace("cyt","co").replace("\n", "")
    return header, gene_identity, start, stop, direction


#### EXTRACT PROTEIN CODING GENES ####
if Find_Genes:
    ###CONFIGS###
    gene_output_dir = file_output_dir + "protein_coding_genes/"
    if not os.path.exists(gene_output_dir):
        os.mkdir(gene_output_dir)
    curren_genes_output_dir = gene_output_dir + "Annelida/"
    if not os.path.exists(curren_genes_output_dir): #creates a output path for all mitochondrial nodes
        os.mkdir(curren_genes_output_dir)
    mitochondrial_genes = ("cox1", "nad3", "nad1", "cox2", "atp8", "cox3", "nad6", "cob", "atp6", "nad5", "nad4l", "nad4", "nad2")
    flanking = 0 #if you would like to grab flanking regions
    if "-f" in options:
        flanking = 21
    if Verbose:
        print (mess_with_font.Green + "\nLocating all protein coding genes using gff annotations\n" + mess_with_font.ENDC, mess_with_font.Blue + str(mitochondrial_genes) + mess_with_font.ENDC)
        if "-f" in options:
            print ("NOTE!!!! --> Grabbing additional 21 nucleotides from flanking regions of genes\n")

    #Obtain gffs and fastas using the extract_file_names definition
    gffs = extract_file_names("gff", GFF_dir)
    fastas = extract_file_names(".fa", FASTA_dir)
    #create a dictionary of fastas and matching annotations:
    fa_gff_dict = {} ####### <-------- The file we will use for alignments and annotations


    for fa in fastas:
        for gff in gffs:
            tmp_gff,tmp_ext = gff.split(".")
            tmp_gff = tmp_gff.replace("_", "")
            tmp_gff = tmp_gff.replace("-", "")
            fatmp = fa.replace("_", "")
            fatmp = fatmp.replace("-", "")
            if tmp_gff in fatmp:
                fa_gff_dict[fa] = gff #fills the dictionary with mito fasta and matching gff

    print (fa_gff_dict)
    for current_fasta, current_gff in fa_gff_dict.items():
        g = open(GFF_dir + current_gff, 'r')
        if Verbose:
            print (mess_with_font.Green + "\nDetails for each gene in the file\n" + mess_with_font.ENDC )
        for line in g.readlines():
            if "trn" in line: #skips if line is not a protein coding gene
                continue
            elif line[0] == "#":
                continue
            elif "CDS" in line:
                continue
            else:
                header, gene_identity, start, stop, direction = parse_gff(line) #grabs relavent details from gff
                if Verbose:
                    print(header, gene_identity, start, stop, direction)
                start = int(start); stop = int(stop)
            for current_gene in mitochondrial_genes: #extracts pcg and writes to appendable file
                if current_gene in gene_identity:
                    gene_file = curren_genes_output_dir + current_gene + ".fasta"
                    extract_region(FASTA_dir + current_fasta, "e", start-flanking, stop+flanking, gene_file, direction, gene_identity + current_fasta)

####Run Muscle on Concatenated File And Each Gene Individual ####
if Alignment:
    #Configs
    align_dir = file_output_dir + "alignments/"
    if not os.path.exists(align_dir):
        os.mkdir(align_dir)

    alignment_output_dir = align_dir + "annelida/"
    if not os.path.exists(alignment_output_dir):
        os.mkdir(alignment_output_dir)


    if Verbose:
        print (mess_with_font.Blue + "\nRunning alignment for each protein coding gene multifasta in", gene_output_dir + mess_with_font.ENDC)

    multi_fastas = extract_file_names(".fasta", curren_genes_output_dir) #extract filenames
    for multi_fasta in multi_fastas:
        multi_fasta = curren_genes_output_dir + multi_fasta
        tmp1, tmp2, tmp3, name = multi_fasta.split("/"); name2, ext = name.split(".") # grab just the name of the gene
        muscle_alignment(multi_fasta,alignment_output_dir + name2 + ".clw") #complete alignment
print ("\n\n{0}Thank{1}{2} You{1} {0}For{1} {2}Using{1} {0}The{1} {2}Program!{1}\n{0}Happy{1} {2}Holidays{1}{0}!!!!{1}\n\n".format (mess_with_font.Green, mess_with_font.ENDC, mess_with_font.Red))
