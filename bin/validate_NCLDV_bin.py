import argparse
import pandas as pd

def validate_NCLDV_bin_from_assessment_summary(ncldv_assessment_summary, hallmark_gene_summary, core_gene_list, bin_id):
    # Load TSV files
    ncldv_df = pd.read_csv(ncldv_assessment_summary, sep="\t")
    hallmark_df = pd.read_csv(hallmark_gene_summary, sep="\t")

    # Calculate the ratio of NCLDV_score = 0
    ncldv_zero_ratio = (ncldv_df["NCLDV_score"] == 0).mean()

    # Count the number of core genes in the bin
    core_gene_count = hallmark_df[core_gene_list].sum().sum()

    # Create the result DataFrame
    result_df = pd.DataFrame({
        "ID": [bin_id],
        "NCLDV_score_zero_ratio": [ncldv_zero_ratio],
        "core_gene_count": [core_gene_count],
        "validate": [0 if ncldv_zero_ratio >= 0.9 and core_gene_count == 0 else 1]
    })
    
    return result_df

def main():
    parser = argparse.ArgumentParser(description="This script validates NCLDV bins using an NCLDV assessment tools summary, a hallmark genes hmmsearch results summary, and a list of NCLDV core genes.")
    parser.add_argument("-b", "--bin_id", required=True, type=str, help="The bin id")
    parser.add_argument('-a', '--ncldv_assessment_summary', required=True, type=str, help="Path to the NCLDV binary assessment summary file in TSV format.")
    parser.add_argument('-m', '--hallmark_gene_summary', required=True, type=str, help="Path to the hallmark gene HMMER results summary file in TSV format.")
    parser.add_argument('-c', '--core_genes', required=True, type=str, help="Comma-separated list of NCLDV core genes.")
    parser.add_argument('-o', '--output', required=True, type=str, help='Path to the output file')

    args = parser.parse_args()

    core_gene_list = args.core_genes.split(',')

    df = validate_NCLDV_bin_from_assessment_summary(args.ncldv_assessment_summary, args.hallmark_gene_summary, core_gene_list, args.bin_id)

    df.to_csv(args.output, index=False, sep='\t')

if __name__ == "__main__":
    main()
