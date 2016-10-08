#!/usr/bin/python                                                                                                                                                                                                                                                       

# """                                                                                                                                                                                                                                                                   
#Usage: converter.py foo.gbk

#Produces a converted file: foo_converted.faa
#"""                                                                                                                                                                                                                                                                    
import sys, os
from Bio import GenBank, SeqIO

gbk_filename = sys.argv[1]
root_name = os.path.splitext(gbk_filename)[0]
faa_filename = root_name + "_converted.faa"

input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print "Dealing with GenBank record %s" % seq_record.id
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            assert len(seq_feature.qualifiers['translation'])==1
            output_handle.write(">%s from %s\n%s\n" % (
                   seq_feature.qualifiers['locus_tag'][0],
                   seq_record.name,
                   seq_feature.qualifiers['translation'][0]))

output_handle.close()
input_handle.close()
print "Done"
