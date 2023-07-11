import argparse
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import iqr
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description='Perform quality control on bin data.')
    parser.add_argument('-f','--fasta_file', required=True, type=str, help='Path to the input FASTA file.')
    parser.add_argument('-c','--coverage_file', required=True, type=str, help='Path to the coverage TSV file.')
    parser.add_argument('-t','--tetranucleotide', required=True, type=str, help='Path to the tetranucleotide frequency TSV file (cgat output).')
    parser.add_argument('-a','--ncldv_assessment_file', required=True, type=str, help='Path to the NCLDV assessment TSV file.')
    parser.add_argument('-o','--output_file', required=True, type=str, help='Path to the output FASTA file.')
    parser.add_argument('-d','--output_data_file', required=True, type=str, help='Path to the output data TSV file.')
    args = parser.parse_args()

    # Load the data
    fasta_sequences = SeqIO.parse(open(args.fasta_file),'fasta')
    coverage_data = pd.read_csv(args.coverage_file, sep='\t')
    tetra_df = pd.read_csv(args.tetranucleotide, sep='\t').T
    new_header = tetra_df.iloc[0]
    tetra_df = tetra_df[1:]
    tetra_df.columns = new_header
    tetra_df.index.name = "ID"
    tetra_df = tetra_df.loc[:, ~(tetra_df == 0).all()]

    ncldv_assessment_data = pd.read_csv(args.ncldv_assessment_file, sep='\t', header=0).set_index('ID')

    # Process coverage data
    coverage_data.set_index('contigName', inplace=True)
    coverage_data.index.name = "ID"
    coverage_data['outlier'] = (coverage_data['totalAvgDepth'] < (coverage_data['totalAvgDepth'].quantile(0.25) - 1.5*iqr(coverage_data['totalAvgDepth']))) | (coverage_data['totalAvgDepth'] > (coverage_data['totalAvgDepth'].quantile(0.75) + 1.5*iqr(coverage_data['totalAvgDepth'])))

    # Process tetranucleotide frequency data
    pca = PCA(n_components=1)
    tetra_df['PC1'] = pca.fit_transform(tetra_df)
    tetra_df['outlier'] = (tetra_df['PC1'] < -2.5*np.std(tetra_df['PC1'])) | (tetra_df['PC1'] > 2.5*np.std(tetra_df['PC1']))

    # Create a new DataFrame for output
    output_data = pd.DataFrame(index=coverage_data.index)
    output_data['TotalAvgDepth'] = coverage_data['totalAvgDepth']
    output_data['depth_outlier'] = coverage_data['outlier'].astype(int)
    output_data['PC1'] = tetra_df['PC1']
    output_data['PC_outlier'] = tetra_df['outlier'].astype(int)
    output_data['NCLDV_score'] = ncldv_assessment_data['NCLDV_score']
    output_data['outlier'] = (output_data[['depth_outlier', 'PC_outlier']].any(axis=1) & (output_data['NCLDV_score'] == 1)).astype(int)

    # Identify outlier contigs
    outlier_contigs = set(output_data[output_data['outlier'] == 1].index)

    # Exclude outlier contigs from the fasta file
    clean_sequences = [fasta for fasta in fasta_sequences if fasta.id not in outlier_contigs]

    # Write the clean sequences to a new fasta file
    SeqIO.write(clean_sequences, args.output_file, "fasta")

    # Write the output data to a new TSV file
    output_data.to_csv(args.output_data_file, sep='\t')

if __name__ == '__main__':
    main()
