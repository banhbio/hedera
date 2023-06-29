"""
non-NCLDV bin definition:
1)	the bin does not have any of the five NCLDV hallmark genes (i.e., MCP, POLB, TFIIS, VLTF3, pATPase);
AND
2)   >90% contigs of the bin were assigned with the GV characteristic score of ‘0’

# input files: 'merged_data_scored.csv' generated from 6.2  'all_putative_bins_contig_list.csv' generated from 7.1 and '8_genes_hits.csv'
"""
import pandas as pd

# Specify the directory where your files are located
file_path = "/Users/yue/Desktop/test/"

# Load 'merged_data_scored.csv'  'all_putative_bins_contig_list.csv' and '8_genes_hits.csv'
score_df = pd.read_csv(file_path + 'merged_data_scored.csv', header=None, skiprows=1)
contig_bin_df = pd.read_csv(file_path + 'all_putative_bins_contig_list.csv', header=None, skiprows=1)
eight_genes_df = pd.read_csv(file_path + '8_genes_hits_modified.csv', header=None, skiprows=1)

# Merge the DataFrames based on the first column
merged_df = pd.merge(score_df, contig_bin_df, on=score_df.columns[0], how="outer")
merged_df = pd.merge(merged_df, eight_genes_df, on=score_df.columns[0], how="outer")

# Rename the columns in the merged DataFrame
merged_df.columns = ["contig", "virsorter2", "CAT", "viralrecall", "149hmm", "GV_score", "bin_id","hallmark_gene_hits"]

## Assess GV density of the bins
# Convert 'GV_score' column to numeric
merged_df['GV_score'] = pd.to_numeric(merged_df['GV_score'], errors='coerce')

# Group merged_df by 'bin_id'
grouped_df = merged_df.groupby('bin_id')

# Calculate occurrence frequency of '0' in 'GV_score' within each group
occurrence_frequency = grouped_df['GV_score'].apply(lambda x: (x == 0).sum())

# Calculate frequency of entries within each group
total_frequency = grouped_df['GV_score'].count()

# Calculate percentage of '0' occurrences within each group
percentage_zero = (occurrence_frequency / total_frequency) * 100

# Create a new DataFrame to store the results
result_df = pd.DataFrame({'bin_id': total_frequency.index, 'Percentage_Zero': percentage_zero.values})

# Merge the result with the original merged_df on 'bin_id'
merged_df = pd.merge(merged_df, result_df, on='bin_id')

# Identify bins meets non-NCLDV criteria
filtered_df = merged_df[(merged_df['Percentage_Zero'] > 90) & (~merged_df['hallmark_gene_hits'].str.contains('MCP|polB|TFIIS|VLTF3|pATPase', na=False))]

# Get unique non-NCLDV bin_id values
unique_bin_ids = filtered_df['bin_id'].unique()

# Write output file
output_file = "/Users/yue/Desktop/test/NonNCLDV_bins_list.txt"
with open(output_file, 'a') as file:
    for bin_id in unique_bin_ids:
        file.write(bin_id + '\n')
