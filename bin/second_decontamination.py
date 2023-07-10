import argparse
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import iqr
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description='Perform quality control on bin data.')
    parser.add_argument('-f','--fasta_file', type=str, help='Path to the input FASTA file.', required=True)
    parser.add_argument('-c','--coverage_file', type=str, help='Path to the coverage TSV file.', required=True)
    parser.add_argument('-t','--tetra_freq_file', type=str, help='Path to the tetranucleotide frequency TSV file (cgat output).', required=True)
    parser.add_argument('-o','--output_file', type=str, help='Path to the output FASTA file.', required=True)
    args = parser.parse_args()

    # Load the data
    fasta_sequences = SeqIO.parse(open(args.fasta_file),'fasta')
    coverage_data = pd.read_csv(args.coverage_file, sep='\t')
    tetra_freq_data = pd.read_csv(args.tetra_freq_file, sep='\t').set_index('kmer').T

    # Process coverage data
    coverage_data['outlier'] = (coverage_data['totalAvgDepth'] < (coverage_data['totalAvgDepth'].quantile(0.25) - 1.5*iqr(coverage_data['totalAvgDepth']))) | (coverage_data['totalAvgDepth'] > (coverage_data['totalAvgDepth'].quantile(0.75) + 1.5*iqr(coverage_data['totalAvgDepth'])))

    # Process tetranucleotide frequency data
    pca = PCA(n_components=1)
    tetra_freq_data['PC1'] = pca.fit_transform(tetra_freq_data)
    tetra_freq_data['outlier'] = (tetra_freq_data['PC1'] < -2.5*np.std(tetra_freq_data['PC1'])) | (tetra_freq_data['PC1'] > 2.5*np.std(tetra_freq_data['PC1']))

    # Identify outlier contigs
    outlier_contigs = set(coverage_data[coverage_data['outlier']]['contigName']).union(set(tetra_freq_data[tetra_freq_data['outlier']].index))

    # Exclude outlier contigs from the fasta file
    clean_sequences = [fasta for fasta in fasta_sequences if fasta.id not in outlier_contigs]

    # Write the clean sequences to a new fasta file
    SeqIO.write(clean_sequences, args.output_file, "fasta")

if __name__ == '__main__':
    main()
