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
genomedir = "/aptmp/ideas2/ivyyue/Ura_1718_42_samples/binning_MetaBAT2/putative_NCLDV_filter/all_individual_assembled_bins_copy"
hmmdir = "/aptmp/ideas2/ivyyue/Ura_1718_42_samples/binning_MetaBAT2/putative_NCLDV_filter/all_individual_bins_hmm_result"

# make a genomelist for output
# calculate genome size of bins
genome_list = {}
file_list=os.listdir(genomedir)
for file_name in file_list:
    if "fa" in file_name:
        print(file_name,end="\r")
        corename = file_name.split(".")[0]+"."+file_name.split(".")[1] #for the case that bin fa files has two dots in filename
        genomesize = 0
        record = SeqIO.to_dict(SeqIO.parse(genomedir+ "/" + file_name, "fasta"))
        for contig in record.keys():
            genomesize += len(record[contig])
        genome_list[corename] = genomesize

# input NCVOG and hmm results
NCVOGs_list = ["NCVOG0022","NCVOG0023","NCVOG0037","NCVOG0038","NCVOG0052","NCVOG0076","NCVOG0236","NCVOG0249","NCVOG0261","NCVOG0262","NCVOG0271","NCVOG0272","NCVOG0273","NCVOG0274","NCVOG0276","NCVOG1060","NCVOG1117","NCVOG1127","NCVOG1164","NCVOG1353"]

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
                        if df[NCVOGs_list[i]][genome] == 0:
                            df[NCVOGs_list[i]][genome] += 1
    i += 1

# weight calculation
weightdict = {"NCVOG0022":0.9,"NCVOG0023":1.1,"NCVOG0037":0.5,"NCVOG0038":1.1,"NCVOG0052":0.9,"NCVOG0076":1,"NCVOG0236":0.8,"NCVOG0249":1,"NCVOG0261":0.7,"NCVOG0262":1,"NCVOG0271":0.9,"NCVOG0272":1,"NCVOG0273":0.8,"NCVOG0274":0.9,"NCVOG0276":0.8,"NCVOG1060":0.6,"NCVOG1117":0.7,"NCVOG1127":0.4,"NCVOG1164":1,"NCVOG1353":0.8}
while i < len(NCVOGs_list):
    df[NCVOGs_list[i]] = df[NCVOGs_list[i]].map(lambda x: x**weightdict[NCVOGs_list[i]])
    i += 1

sumcol = df.sum(axis=1)
df["weight"] = sumcol
df['genomesize'] = df.index.map(genome_list)
df["final_weight"] = df["weight"]/(np.log10(df['genomesize'])-4)

# output a NCVOG distribution and weight file
df.to_csv("NCVOG_genome_binary_distribution_weight.tsv", sep='\t')