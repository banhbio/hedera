#!/usr/bin/env/ python
import argparse
import pandas as pd
import numpy as np


def process_hmm_results(hmm_file, hallmarks_list, scoredict):
    df = pd.DataFrame(np.zeros((1, len(hallmarks_list)), dtype=int), columns=hallmarks_list)

    with open(hmm_file, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            if line.startswith("#"):
                continue
            else:
                data = line.split()
                hallmark = data[2]
                score = float(data[5])
                if score >= scoredict[hallmark]:
                    df[hallmark] += 1
    
    df['total'] = df.sum(axis=1)
    return df

def main():
    parser = argparse.ArgumentParser(description='Detect NVLDVs hallmark genes from hmmsearch results.')
    parser.add_argument('-t', '--hmm', type=str, help='Path to the hmmsearch results file (tblout)')
    parser.add_argument('-o', '--output', type=str, help='Path to the output file')
    parser.add_argument('-n', '--nhallmark', type=str, help='Comma-separated list of NCLDV hallmark genes')
    parser.add_argument('-s', '--score', type=str, help='Comma-separated list of score thresholds for each NCLDV hallmark gene')

    args = parser.parse_args()

    hallmarks_list = args.nhallmark.split(',')
    score_list = [float(i) for i in args.score.split(',')]
    scoredict = dict(zip(hallmarks_list, score_list))

    df = process_hmm_results(args.hmm, hallmarks_list, scoredict)

    df.to_csv(args.output, index=False, sep='\t')

if __name__ == "__main__":
    main()
