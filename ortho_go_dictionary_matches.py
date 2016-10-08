import pickle, sys, re

pickle_loc = '/home/genome/joseph7e/functional_16S_stuff/go_dict.pickle'
pickle_dict2 = ''

accession_dict = {} #accession:[match,unmatch]
print ("building accession dictionary")


accession_file = open('/home/genome/joseph7e/gene_16S/analysis_orthofinder/cluster_99_accession_matches', 'r')
for line in accession_file:
    line = line.rstrip()
    accession, match, unmatch = line.split("\t")
    accession_dict[accession] = [match,unmatch]



bio_dict = {} #P --> go: [accessions]
mol_dict = {} #F --> go: [accessions]
cell_dict = {} #C --> go: [accessions]


print ("Importing database dictionary with pickle")

with open(pickle_loc, 'rb') as handle:
   dict1 = pickle.load(handle)


for id, lists in dict1.items():
    accessions = lists[0]
    gos = lists[1]
    if len(gos) > 0:
        bios = re.findall(r"P:(.*?)[,']", str(gos))
        mols = re.findall(r"F:(.*?)[,']", str(gos))
        cells = re.findall(r"C:(.*?)[,']", str(gos))
        for b in bios:
            if b in bio_dict.keys():
                bio_dict[b].extend(accessions)
            else:
                bio_dict[b] = accessions
        for m in mols:
            if m in mol_dict.keys():
                mol_dict[m].extend(accessions)
            else:
                mol_dict[m] = accessions
        for c in cells:
            if c in cell_dict.keys():
                cell_dict[c].extend(accessions)
            else:
                cell_dict[c] = accessions

final_bio_dict = {}
final_mol_dict = {}
final_cell_dict = {}

for acc, score in accession_dict.items():
    for go, accs in bio_dict.items():
        for a in accs:
            if acc in a:
                if go in final_bio_dict.keys():
                    final_bio_dict[go][0] += int(score[0])
                    final_bio_dict[go][1] += int(score[1])
                else:
                    final_bio_dict[go] = [int(score[0]), int(score[1])]
    for go, accs in cell_dict.items():
        for a in accs:
            if acc in a:
                if go in final_cell_dict.keys():
                    final_cell_dict[go][0] += int(score[0])
                    final_cell_dict[go][1] += int(score[1])
                else:
                    final_cell_dict[go] = [int(score[0]), int(score[1])]
    for go, accs in mol_dict.items():
        for a in accs:
            if acc in a:
                if go in final_mol_dict.keys():
                    final_mol_dict[go][0] += int(score[0])
                    final_mol_dict[go][1] += int(score[1])
                else:
                    final_mol_dict[go] = [int(score[0]), int(score[1])]

output_bio = open('gene_ontology_enrichment_Bio.tsv','w')
output_mol = open('gene_ontology_enrichment_Mol.tsv','w')
output_cell = open('gene_ontology_enrichment_Cell.tsv','w')


def save_dict_stats(dictionary, out):
    for key, value in dictionary.items():
        out.writelines(key + '\t' + str(value[0]) + '\t' + str(value[1]))

save_dict_stats(final_bio_dict, output_bio)
save_dict_stats(final_mol_dict, output_mol)
save_dict_stats(final_cell_dict, output_cell)


#awk -F"\t" '{print $2 "\t" $3}' gene_ontology_enrichment_Bio.tsv | sort -nk 2