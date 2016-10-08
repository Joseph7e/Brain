import re, sys
from collections import defaultdict


def getReference(f):
    return "".join(f.readlines()[1:]).replace("\n", "")

def parse(file, ref):
    line = file.readline()
    haplo = defaultdict(dict)
    haplo_index = {}
    while line[0:2] == "##":
        line = file.readline()

    counter = 0
    for sample in line.split("\t")[9:]:
        sample = re.search("rev(\d*)", sample).group(1)
        haplo_index[counter] = sample

        for i in range(4):
            haplo[sample][i] = bytearray(ref)

        counter += 1

    line = file.readline()
    while line:
        elements = line.split('\t')
        pos = int(elements[1]) - 1
        orig = elements[3]
        snp = elements[4]

        for index, h in enumerate(elements[9:]):
            hp = haplo[haplo_index[index]]
            for allele_num, allele in enumerate(h.split(":")[0].split("/")):
                if int(allele) == 0:
                    continue

                if len(orig) > 1:
                    hp[allele_num][pos + 1] = '-'
                else:
                    hp[allele_num][pos] = snp

                # for snp_pos, s in enumerate(snp):
                #     hp[allele_num][pos + snp_pos] = s
                #     hp[allele_num][pos + snp_pos] = s

        # print orig, snp, ref[int(pos) - 1]

        line = file.readline()
    stuff = set()
    for hap, alleles in haplo.iteritems():
        temp_set = set()
        succ_count = 0
        out = open("rev_" + str(hap) + ".fa", "w")


        for allele_num, allele in alleles.iteritems():
            stuff.add(str(allele))
            success = (str(allele)) not in temp_set
            temp_set.add(str(allele))
            if success:
                succ_count += 1
                out.write(">" + "rev_" + str(hap) + "_" + str(succ_count) + "\n" + allele + "\n")
        out.close()

    out = open("unique_alleles.fa", "w")
    print "\nUNIQUE ALLELES\n"
    for idx,i in enumerate(stuff):
        out.write(">unique_" + str(idx + 1) + "\n" + i + "\n")
    out.close()

if __name__ == '__main__':
    ref = None
    with open(sys.argv[1], "r") as f:
        ref = getReference(f)
    with open(sys.argv[2], "r") as f:
        parse(f, ref)
