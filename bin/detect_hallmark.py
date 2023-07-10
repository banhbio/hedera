#!/usr/bin/env/ python
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO


def process_hmm_results(contig_list, hmm_file, hallmarks_list, scoredict):
    # Initialize an empty DataFrame for total counts
    df_total = pd.DataFrame(np.zeros((1, len(hallmarks_list)), dtype=int), columns=hallmarks_list)
    # Initialize an empty DataFrame for individual contig counts
    df_contigs = pd.DataFrame(np.zeros((len(contig_list), len(hallmarks_list)), dtype=int), columns=hallmarks_list, index=contig_list)

    with open(hmm_file, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            if line.startswith("#"):
                continue
            else:
                data = line.split()
                pep = data[0]
                contig = "_".join(pep.split("_")[0:-1])
                hallmark = data[2]
                score = float(data[5])
                if score >= scoredict[hallmark]:
                    df_total[hallmark] += 1
                    df_contigs.loc[contig, hallmark] += 1

    df_total['total'] = df_total.sum(axis=1)
    df_contigs['total'] = df_contigs.sum(axis=1)
    return df_total, df_contigs

def main():
    parser = argparse.ArgumentParser(description='Detect NVLDVs hallmark genes from hmmsearch results.')
    parser.add_argument('-f', '--fasta', type=str, help='Path to the fasta file')
    parser.add_argument('-t', '--hmm', type=str, help='Path to the hmmsearch results file (tblout)')
    parser.add_argument('-o', '--total_output', type=str, help='Path to the output (total count) TSV file')
    parser.add_argument('-O', '--contigs_output', type=str, help='Path to the output (count per contigs) TSV file')
    parser.add_argument('-n', '--nhallmark', type=str, help='Comma-separated list of NCLDV hallmark genes')
    parser.add_argument('-s', '--score', type=str, help='Comma-separated list of score thresholds for each NCLDV hallmark gene')

    args = parser.parse_args()

    contig_list = [record.id for record in SeqIO.parse(args.fasta, "fasta")]
    hallmarks_list = args.nhallmark.split(',')
    score_list = [float(i) for i in args.score.split(',')]
    scoredict = dict(zip(hallmarks_list, score_list))

    df_total, df_contigs = process_hmm_results(contig_list, args.hmm, hallmarks_list, scoredict)

    df_total.to_csv(args.total_output, index=False, sep='\t')
    df_contigs.to_csv(args.contigs_output, index=True, index_label="ID", sep='\t')

if __name__ == "__main__":
    main()
