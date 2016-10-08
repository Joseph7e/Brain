import json
from pprint import pprint
import sys
my_json = sys.argv[1]#'/home/genome/kdiaz/ALL_CONTIGS/sample2/BlobDB.json'

with open(my_json, encoding='utf-8') as data_file:
    data = json.loads(data_file.read())


#pprint(data)
#pprint (data["dict_of_blobs"]["NODE_10465_length_6188_cov_8.55981_ID_4083273"])

for key, value in data["dict_of_blobs"].items():
    print (key,'\t',data["dict_of_blobs"][key]['covs']['cov0'],'\t',data["dict_of_blobs"][key]['gc'],'\t',data["dict_of_blobs"][key]['length'],'\t', data["dict_of_blobs"][key]["taxonomy"]["bestsum"]["superkingdom"]["tax"],'\t', data["dict_of_blobs"][key]["taxonomy"]["bestsum"]["phylum"]["tax"],'\t', data["dict_of_blobs"][key]["taxonomy"]["bestsum"]["order"]["tax"],'\t', data["dict_of_blobs"][key]["taxonomy"]["bestsum"]["family"]["tax"],'\t', data["dict_of_blobs"][key]["taxonomy"]["bestsum"]["genus"]["tax"],'\t', data["dict_of_blobs"][key]["taxonomy"]["bestsum"]["species"]["tax"])

#sed 's/NODE_//' contig_tax_info_sample_a1.tsv | sort -n | sed 's/^/NODE_/' > contig_tax_info_sample_a1.tsv
