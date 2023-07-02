#!/usr/bin/env/ python
import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np

def calculate_genome_size(fasta_file):
    genome_size = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        genome_size += len(record.seq)
    return genome_size

def process_hmm_results(hmm_file, genome_size, NCVOGs_list, weightdict):
    df = pd.DataFrame(np.zeros((1, len(NCVOGs_list))), columns=NCVOGs_list)

#Multiple NCVOGs in one genome count as one.
    with open(hmm_file, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            if line.startswith("#"):
                continue
            else:
                for ncvog in NCVOGs_list:
                    if ncvog in line:
                        df[ncvog] = 1

    for ncvog in NCVOGs_list:
        df[ncvog] = df[ncvog]*weightdict[ncvog]

    df["weight"] = df.sum(axis=1)
    df['genomesize'] = genome_size
    df["final_score"] = df["weight"] / (np.log10(df['genomesize']) - 4)

    return df

def main():
    parser = argparse.ArgumentParser(description='Calculate NVLDVs core genes index from hmmsearch results and genome size.')
    parser.add_argument('-f', '--fasta', type=str, help='Path to the fasta file')
    parser.add_argument('-t', '--hmm', type=str, help='Path to the hmmsearch results file (tblout)')
    parser.add_argument('-o', '--output', type=str, help='Path to the output file')
    parser.add_argument('-n', '--ncvogs', type=str, help='Comma-separated list of NCVOGs')
    parser.add_argument('-w', '--weights', type=str, help='Comma-separated list of weights for each NCVOG')

    args = parser.parse_args()

    NCVOGs_list = args.ncvogs.split(',')
    weight_list = [float(i) for i in args.weights.split(',')]
    weightdict = dict(zip(NCVOGs_list, weight_list))

    genome_size = calculate_genome_size(args.fasta)
    df = process_hmm_results(args.hmm, genome_size, NCVOGs_list, weightdict)

    df.to_csv(args.output, index=False, sep='\t')

if __name__ == "__main__":
    main()
