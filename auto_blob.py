#!/usr/bin/python3
#AUTHOR: Joseph Sevigny and Krystalle Diaz
#Purpose: Run bowtie and blob tools.

import sys, subprocess, os


### Global Varaiables

input_contigs = sys.argv[1] #contigs.fasta
sample_name = sys.argv[2] # Nem-A02
forward_reads = sys.argv[3]
reverse_reads = sys.argv[4]
unpaired_reads = sys.argv[5]

blast_db = '/home/genome/kdiaz/databases/nt'

bash_script = open('bash_run_bowtie_and_blast.sh', 'w')
bowtie_command = "bowtie2 -x "+sample_name+"build"+" -p 30 -1 "+forward_reads+" -2 "+reverse_reads+" -U "+unpaired_reads+" &"
blast_command = "blastn -db /home/genome/kdiaz/databases/nt -query " + input_contigs +  " -outfmt '6 qseqid staxids qcovs' -evalue 1e-10 -num_threads 16 -out tax_file &"

# #index with bowtie
index_command = ['bowtie2-build', input_contigs, sample_name+'build']
si = subprocess.Popen(index_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
si.communicate()

bash_script.writelines("#!/bin/bash\n"+bowtie_command+"\n"+blast_command+"\n"+"wait")
bash_script.close()
os.system("chmod 755 bash_run_bowtie_and_blast.sh")
os.system("./bash_run_bowtie_and_blast.sh")

#align with bowtie --> dont wait
#bowtie_command = ['bowtie2', '-x', sample_name+'build', '-p', '30', '-1', forward_reads, '-2', reverse_reads, '-U', unpaired_reads]
#sb = subprocess.Popen(bowtie_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
samtools_command = ['samtools']


#create blast file for blobtools wait
#os.system("blastn -db /home/genome/kdiaz/databases/nt -query " + input_contigs +  " -outfmt '6 qseqid staxids qcovs' -evalue 1e-10 -num_threads 16 -out tax_file")
# blast_command = ['blastn', '-db', blast_db, '-query', input_contigs, '-outfmt', '"6', 'qseqid', 'staxids', 'qcovs\"', '-evalue', '1e-10', '-num_threads', '16', '-out', 'tax_file']
# sblast = subprocess.Popen(blast_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


#####Blobtools