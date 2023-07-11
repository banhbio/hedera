import argparse
import pandas as pd
import numpy as np

def get_coverages(depth):
    df = pd.read_csv(depth, sep='\t')
    coverages = df['totalAvgDepth'].values
    return coverages

def calculate_coef_of_variance(coverages):
    coef_of_variance = np.std(coverages) / np.mean(coverages)
    return coef_of_variance

def identify_potential_chimera_bins(hallmark_gene_summary, depth, scgs_list, bin_id):
    coverages = get_coverages(depth)
    coef_of_variance = calculate_coef_of_variance(coverages)

    hallmark_df = pd.read_csv(hallmark_gene_summary, sep='\t')
    multiple_scgs = ','.join([gene for gene, count in hallmark_df[scgs_list].sum().to_dict().items() if count > 1])
    delinage_candidate = 1 if coef_of_variance > 0.01 and multiple_scgs else 0

    # Create a new DataFrame
    new_df = pd.DataFrame({
        "ID": [bin_id],
        "coefficient_of_variance": [coef_of_variance],
        "multiple_scgs": [multiple_scgs],
        "delinage": [delinage_candidate]
    })
    
    return new_df

def main():
    parser = argparse.ArgumentParser(description='Process a hallmark gene summary file and a depth file to identify potential chimera bins.')
    parser.add_argument("-b", "--bin_id", required=True, type=str, help="The bin id")
    parser.add_argument('-m', '--hallmark_gene_summary', required=True, type=str, help="Path to the hallmark gene HMMER results summary per contig file in TSV format.")
    parser.add_argument('-s', '--scgs', required=True, type=str, help='Comma-separated list of NCLDV single copy genes.')
    parser.add_argument('-d', '--depth', required=True, type=str, help="Path to the depth file (totalAvgDepth in 3rd column) in TSV format.")
    parser.add_argument('-o', '--output', required=True,type=str, help='Path to the output result TSV file.')
    args = parser.parse_args()
    
    scgs_list = args.scgs.split(',')

    df = identify_potential_chimera_bins(args.hallmark_gene_summary, args.depth ,scgs_list, args.bin_id)
    
    df.to_csv(args.output, index=False, sep='\t')

if __name__ == "__main__":
    main()