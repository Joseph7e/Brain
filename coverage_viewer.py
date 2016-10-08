import numpy as np
import sys


if __name__ == '__main__':
    mito_filename = sys.argv[1]
    sam_filename = sys.argv[2]
    header = ""
    contigs = {}
    base_arrays = {}
    sequence = bytearray()

    # read contents of mitochondrial file and store in contig dictionary
    with open(mito_filename, "r") as f:
        block = f.readlines(50000)
        seq_length = 0
        while block:
            for line in block:
                if line[0] == ">":
                    if header:
                        base_arrays[header] = np.zeros(seq_length, dtype=int)
                        contigs[header] = sequence
                    header = line.rstrip()[1:]
                    seq_length = 0
                    sequence = bytearray()
                else:
                    seq_length += len(line) - 1
                    sequence.extend(line.rstrip())


            block = f.readlines(50000)

        contigs[header] = sequence
        base_arrays[header] = np.zeros(seq_length + 1, dtype=int)



    with open(sam_filename, "r") as f:
        block = f.readlines(50000)
        seq_length = 0
        while block:
            for line in block:
                elements = line.split("\t")
                # don't parse unmapped reads
                if elements[2] == "*":
                    continue

                start = int(elements[3])
                contig = elements[2]
                rev_comp = int(elements[1]) & 16
                if not rev_comp:
                    base_arrays[contig][start:start + 250] += 1
                else:
                    base_arrays[contig][start - 250:start] += 1

            block = f.readlines(50000)



    threshold = 10
    for key in contigs:
        contig_sequence = contigs[key]
        base_array = base_arrays[key]


        print(contig_sequence)
        for i in range(np.amax(base_array) / threshold + 1):
            i = (i * threshold) + 1
            print "".join(map(lambda x: "*" if x >= i else " ", base_array))
