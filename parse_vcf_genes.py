import sys
import re
import argparse
from collections import namedtuple, defaultdict

# record containing some of the columns found in a GFF file
# a named tuple makes it easier to remember what's stored in it
Feature = namedtuple('Feature', 'contig, type, start, stop, direction, notes')


# overloaded error method so it prints help info when user gives wrong arguments
class NewParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('ERROR: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


def create_parser():
    """
    Creates a help parser (using custom argparser class)
    :return:
    """
    p = NewParser()

    p.add_argument('reference', type=str,
                   help = "Fasta reference file that reads were mapped to.")

    p.add_argument('gff', type=str,
                   help = "GFF file containing reference genome annotations.")

    p.add_argument('vcf', type=str,
                   help = "VCF file to parse.")

    args = p.parse_args(sys.argv[1:])
    return args


def create_reference(filename):
    """
    Quick and dirty creation of contigs dictionary from file.
    """
    contigs = {}
    with open(filename, "r") as f:
        header = ""
        seq = ""
        for line in f:
            if line[0] == ">":
                if header:
                    contigs[header] = seq
                header = line[1:].rstrip()
                seq = ""
            else:
                seq += line.rstrip()

    return contigs

def parse_gff(filename):
    """
    :param filename: location to GFF file to open
    :return:
        A dictionary where the keys are the names of the contigs, and the
        values are the GFF annotations found on the contig.
    """
    gff_contigs = defaultdict(list)
    with open(filename, "r") as f:
        for line in f:
            if line == "##FASTA\n": # filter out fasta crap at bottom of gff
                break
            if line[0] == "#": # filter out comments from gff
                continue


            elements = line.rstrip().split("\t")

            # this record is used to store GFF column info and append it to list
            entry = Feature(elements[0], elements[2], elements[3], elements[4], elements[6], elements[8])
            gff_contigs[entry.contig].append(entry)

    return gff_contigs



def extract_feature_notes(contig, type, start, stop, direction, notes, pos, ref, alt):
    """
    Prints information found within the notes of a feature (if available).
    """
    notes = re.sub("%\d+", " ", notes)
    desc = re.search("Note=(.*?);", notes)
    id = re.search("ID=(.*?);", notes)
    name = re.search("Name=(.*?);", notes)

    if id and name and desc:
        print ('Contig: ' + contig)
        print ('Annotation_type: ' + type)
        print ('Start:Stop:Direction '+ start+':'+stop+':'+direction)
        print ("Position: " + pos)
        print ('Ref:Alt: ' + ref + ':'+ alt)
        print("ID: " + id.group(1))
        # print("Name: " + name.group(1))
        print("Description: " + desc.group(1))
    # if name:
    #     print("Name: " + name.group(1))
    #
    # if desc:
    #     print("Desc: " + desc.group(1))


    print()




def parse_vcf(filename, gff_contigs):
    """
    Iterates over a vcf file.

    :param filename:  name of the VCF file to parse
    :param gff_contigs: output of parse_gff. Used to find snps in genes.
    """
    with open(filename, "r") as f:
        for line in f:
            if line[0] == "#": # ignore comment lines in vcf file
                continue

            elements = line.rstrip().split("\t")
            name = elements[0] # name of the contig / chromosome SNP found on
            pos = elements[1] # position of SNP in contig
            ref = elements[3] # reference basepair(s)
            alt = elements[4] # SNP basepair(s) at same location
            qual = float(elements[5]) # SNP quality, might want to filter this
            contig = gff_contigs[name] # grab contig where SNP is located

            # iterate over GFF annotations on his contig, printing those that
            # overlap with the position of the SNP.
            for feature in contig:
                if pos >= feature.start and pos <= feature.stop:
                    extract_feature_notes(feature.contig, feature.type, feature.start, feature.stop, feature.direction, feature.notes, pos, ref, alt) # prints info


if __name__ == '__main__':

    # creates help parser (prints usage info on error or -h)
    args = create_parser()

    # parses fasta file into contigs dictionary
    contigs = create_reference(args.reference)

    # parses GFF file into dictionary of contigs, value is list of annotations
    gff_contigs = parse_gff(args.gff)

    # iterate over VCF entries (SNPs) and do stuff using gff_contigs
    parse_vcf(args.vcf, gff_contigs)


