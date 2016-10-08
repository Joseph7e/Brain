import os, sys, gzip, re

### R for GO tern analysis , Bioconductor

#

#PURPOSE: inputs a file of accessions, one per line, outputs file with go terms added
#you will need to run twice with both databases to get all

accession_file = open(sys.argv[1], 'r')
database = gzip.open(sys.argv[2], 'rt') #gzipped database file
output = open(sys.argv[3], 'w')

print ('current database = ', sys.argv[2])

accession_dict = {} #accession:[bio,mol,cell]
print ("building accession dictionary from file", sys.argv[1])

for line in accession_file:
    accession = line.rstrip()
    accession_dict[accession] = []
print ('total number of accessions = ', len(accession_dict))

#### Construct Dictionaries
print ("importing data and adding gene ontology information to accession dict")

ID = ''; uniprot = ''; accessions = []; gos = []; cogs = []; kos = []
count = 0; real_count = 0; acc_count = 0

for line in database:
    count += 1; real_count += 1
    if count == 50000:
        print ('current line in database', real_count, 'and', acc_count, "accessions found")
        count = 0
    if line.startswith("ID"):
        if len(accessions) > 0:# and len(gos) > 0:
            # for acc, score in accession_dict.items():
            for a in accessions:
                if a in accession_dict.keys():
                    if accession_dict[a] == [[],[],[],[]] or accession_dict[a] == []:
                        acc_count += 1
                        bios = []; mols = []; cells = []
                        for g in gos:
                            if g.startswith('P:'):
                                bios.append(g[2:])
                            if g.startswith('F:'):
                                mols.append(g[2:])
                            if g.startswith('C:'):
                                cells.append(g[2:])
                        # bios = re.findall(r"P:(.*?)[;]", str(gos))
                        # mols = re.findall(r"F:(.*?)[;]", str(gos))
                        # cells = re.findall(r"C:(.*?)[;]", str(gos))
                        accession_dict[a] = [bios, mols, cells, cogs]
                        #print (a, accession_dict[a])
                        output.write(a + '\t' + str(bios).lstrip().rstrip().replace(" ", "_") + '\t' + str(mols).lstrip().rstrip().replace(" ", "_") + '\t' + str(cells).lstrip().rstrip().replace(" ", "_") + '\t' + str(cogs) + '\t' + str(kos) + '\n')
                    else:
                        print ('got one that was repeated', accessions)
        accessions = []; gos = []; cogs = []; kos = []
        ID = re.findall(r"ID\s*(.*?)\s", line)[0]
    if line.startswith("DR") and "RefSeq" in line:
        tmp, acc_temp, tmp = line.split(';')
        accession = acc_temp.lstrip()
        if '.' in accession:
            accessions.append(accession[:-2])
        else:
            accessions.append(accession)
        # refseq = re.findall(r"RefSeq;\s*(.*?)[;.]", line)
        # for ttt in refseq:
        #     if '.' in ttt:
        #         accessions.append(ttt[:-2])
        #     else:
        #         accessions.append(ttt)
    if line.startswith("DR") and "GO;" in line:
        go = re.findall(r"([PCF]:.*);", line)
        try:
            gos.append(go[0])
        except IndexError:
            continue
    if line.startswith("DR") and "COG" in line:
        if "eggNOG" in line:
            cog = re.findall(r"(COG....);", line)
            try:
                cogs.append(cog[0])
            except IndexError:
                continue
    if line.startswith("DR") and "KO;" in line:
        ko = re.findall(r"(K0..*?);", line)
        try:
            kos.append(ko[0])
        except IndexError:
            continue


print ("{} of the total {} accessions found in this database".format(str(acc_count), str(len(accession_dict))))
# print ("Saving information to output")
#
# output_bio = open('gene_ontology_enrichment_Bio.tsv','w')
# output_mol = open('gene_ontology_enrichment_Mol.tsv','w')
# output_cell = open('gene_ontology_enrichment_Cell.tsv','w')
#
#
# def save_dict_stats(dictionary, out):
#     for key, value in dictionary.items():
#         out.writelines(key + '\t' + str(value[0]) + '\t' + str(value[1])+ '\n')
#
# save_dict_stats(bio_dict, output_bio)
# save_dict_stats(mol_dict, output_mol)
# save_dict_stats(cell_dict, output_cell)