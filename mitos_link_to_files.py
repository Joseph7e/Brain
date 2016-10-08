#!/usr/bin/python3

#Author: Joseph Sevigny
#Purpose: Provided a link to a mitos output website, download and name all datafiles from link
#Usage: mitos_link_to_files.py <web_url> <out_name>

import os, sys, re

web_url = sys.argv[1]#'http://mitos.bioinf.uni-leipzig.de/result.py?hash=ZUOcm12Y'#sys.argv[1]
output_name = sys.argv[2]#'sample_8'

html = ''
import urllib.request
with urllib.request.urlopen(web_url) as response:
    html = response.read()
    html = str(html)

gff_link = 'http://mitos.bioinf.uni-leipzig.de/' + re.findall(r' href="(download\.py\?type=gff.*)?">GFF file', html)[0]
faa_link = 'http://mitos.bioinf.uni-leipzig.de/' + re.findall(r' href="(download\.py\?type=faa.*)?">FAA file', html)[0]
fas_link = 'http://mitos.bioinf.uni-leipzig.de/' + re.findall(r' href="(download\.py\?type=fas.*)?">FAS file', html)[0]

os.system('wget "'+gff_link+'" -O '+output_name+'.gff')
os.system('wget "'+fas_link+'" -O '+output_name+'.fas')
os.system('wget "'+fas_link+'" -O '+output_name+'.fna')
#lines = html.split('\n')



#for line in lines:
#    print (line)
# for line in html.split('\n'):#readlines():
#     print (line)


class Mitos_files(object):
    def __init__(self, html_results):
        for line in html_results.readlines():
            print (line)
        self.sample = ''
        self.gff = ''
        self.faa = ''
        self.fas = ''
