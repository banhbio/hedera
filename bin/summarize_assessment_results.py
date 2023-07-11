#!/usr/bin/env/ python
import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np

def process_results(contigid_list, viralrecall, virsorter2, cat, ncldvogs):
    # Create DataFrame with specified column names
    df = pd.DataFrame(columns=["ID","viralrecall","virsorter2","cat","149ncldvogs"])

    # Set 'ID' column as contigid_list and initialize all other columns with 0
    df['ID'] = contigid_list
    df['viralrecall'] = 0
    df['virsorter2'] = 0
    df['cat'] = 0
    df['149ncldvogs'] = 0

    # Process the viralrecall result file
    with open(viralrecall, 'r') as vr_file:
        next(vr_file)  # Skip header
        for line in vr_file:
            fields = line.strip().split('\t')
            contig_id = fields[1]
            score = float(fields[3])
            if score > 0:
                df.loc[df['ID'] == contig_id, 'viralrecall'] = 1

    # Process the virsorter2 result file
    with open(virsorter2, 'r') as vs_file:
        next(vs_file)  # Skip header
        for line in vs_file:
            fields = line.strip().split('\t')
            contig_id = fields[0].split('||')[0]  # Remove unwanted string after '||'
            classification = fields[7]
            if classification == 'NCLDV':
                df.loc[df['ID'] == contig_id, 'virsorter2'] = 1

    # Process the CAT result file
    with open(cat, 'r') as cat_file:
        next(cat_file)  # Skip header
        for line in cat_file:
            fields = line.strip().split('\t')
            if len(fields) >= 10:  # Ensure the 10th field exists
                contig_id = fields[0]
                value_field = fields[9]
                A, BC = value_field.split(' (') # A (B): C
                B, C = BC.split('): ')
                if A == 'Nucleocytoviricota' and float(C) > 0.50:
                    df.loc[df['ID'] == contig_id, 'cat'] = 1

    # Process the ncldvogs result file
    with open(ncldvogs, 'r') as n_file:
        for line in n_file:
            if not line.startswith('#'):  # Skip comments
                fields = line.strip().split()
                protein_id = fields[0]
                contig_id = '_'.join(protein_id.split('_')[:-1])  # Remove last field after '_' from prodigal protein ID
                if contig_id in df['ID'].values:
                    df.loc[df['ID'] == contig_id, '149ncldvogs'] = 1

    # Add a new column 'NCLDV_score' which is the sum of the values of the columns 2-4
    df['NCLDV_score'] = df.iloc[:, 1:5].sum(axis=1)
    return df

def main():
    parser = argparse.ArgumentParser(description='Calculate NVLDV scores for each contig in single bin.')
    parser.add_argument('-f', '--fasta', required=True, type=str, help='Path to the fasta file')
    parser.add_argument('-r', '--viralrecall', required=True, type=str, help='Path to the vuralrecall result file')
    parser.add_argument('-s', '--virsorter2', required=True, type=str, help='Path to the virsorter2 result file')
    parser.add_argument('-c', '--cat', required=True, type=str, help='Path to the CAT result file')
    parser.add_argument('-n', '--ncldvogs', required=True, type=str, help='Path to the 149 NCLDVOGs hmmsearch results file (tblout)')
    parser.add_argument('-o', '--output', required=True, type=str, help='Path to the output file')

    args = parser.parse_args()

    contigid_list = [record.id for record in SeqIO.parse(args.fasta, "fasta")]
    df = process_results(contigid_list, args.viralrecall, args.virsorter2, args.cat, args.ncldvogs)

    df.to_csv(args.output, index=False, sep='\t')

if __name__ == "__main__":
    main()
