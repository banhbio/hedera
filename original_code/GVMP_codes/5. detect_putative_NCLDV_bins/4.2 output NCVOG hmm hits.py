#contig id in hmm.out need to have bin id in the first place 

#!/usr/bin/env python
# coding: utf-8
import sys
import os
import pandas as pd
import argparse
import numpy as np
from Bio import SeqIO
'''
It's ugly but useful
'''


# fna files of bins
# hmm results
genomedir = "/aptmp/ideas2/ivyyue/tara_ocean_polar_protists_PRJEB9691/PRJEB9691/28_0.8_samples/coassembly/binning/putative_NCLDV_bins_cp"
hmmdir = "/aptmp/ideas2/ivyyue/tara_ocean_polar_protists_PRJEB9691/PRJEB9691/28_0.8_samples/coassembly/binning/putative_NCLDV_filter/NCLDV_bins_hmmsearch"

# make a genomelist for output
# calculate genome size of bins
genome_list = {}
g=os.walk(genomedir)
for path,dir_list,file_list in g:
    for file_name in file_list:
        corename = file_name.split(".")[0]+"."+file_name.split(".")[1] #for the case that bin fa files has two dots in filename
        genomesize = 0
        record = SeqIO.to_dict(SeqIO.parse(path + "/" + file_name, "fasta"))
        for contig in record.keys():
            genomesize += len(record[contig])
        genome_list[corename] = genomesize



# input NCVOG and hmm results
NCVOGs_list = ["DNAPOLB","MCP","Primase","RNA_a","RNA_b","TFIIS","VLTF3","pATPase"]

zeromatrix = np.zeros((len(genome_list.keys()),len(NCVOGs_list)))
df = pd.DataFrame(zeromatrix, index = genome_list.keys(),columns = NCVOGs_list)

i = 0
while i < len(NCVOGs_list):
    with open(f"{hmmdir}/{NCVOGs_list[i]}.out", "r") as infile:
        lines = infile.readlines()
        for line in lines:
            if line.startswith("#"):
                continue
            else:
                for genome in genome_list.keys():
                    if genome + "_" in line:
                            df[NCVOGs_list[i]][genome] += 1
    i += 1
df.to_csv("NCVOG_hmm_hits_in_bins.tsv", sep='\t')