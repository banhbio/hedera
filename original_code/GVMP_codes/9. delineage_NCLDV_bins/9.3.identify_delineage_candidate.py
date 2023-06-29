# Input file: 'hallmark_hmm_hits_in_bins.tsv'
# Identify NCLDV bins with multiple copy of SCG (superfamily II helicase (SFII/TFIIS), viral late transcription factor 3 (VLTF3), A32-ATPase (A32), and B-family DNA Polymerase (PolB)) (Aylward et al., 2021))
import pandas as pd

# Read the file into a DataFrame with stripped whitespace in column names
df = pd.read_csv('/Users/yue/Desktop/test/hallmark_hmm_hits_in_bins.tsv', sep='\s+', index_col=0, header=0, skipinitialspace=True)

# Filter the rows based on the conditions
filtered_df = df[(df['TFIIS'] > 1) | (df['VLTF3'] > 1) | (df['pATPase'] > 1) | (df['DNAPOLB'] > 1)]

# Save the desired line as a CSV file
filtered_df.to_csv('/Users/yue/Documents/masters_project/P002_MAGGING_PIPELINE/MAGGING_PIPELINE_codes/9. delineage_NCLDV_bins/test/multiple_SCG_bins.csv', sep=',')

# Retrieve the NCLDV bin with CV greater than 1

# Read the 'CV_per_bin.csv' generated from 9.1.2.
df2 = pd.read_csv('/Users/yue/Documents/masters_project/P002_MAGGING_PIPELINE/MAGGING_PIPELINE_codes/9. delineage_NCLDV_bins/test/CV_per_bin.csv')

# Filter the rows based on the condition
filtered_df2 = df2[df2.iloc[:, 1] > 1]

# Retrieve the bin with CV greater than 1 AND multiple SCG
common_rows = filtered_df[filtered_df.index.isin(filtered_df2['BIN_ID'])]
common_strings = common_rows.index.tolist()
# Output the candiate bin id list for delineage 
with open('/Users/yue/Documents/masters_project/P002_MAGGING_PIPELINE/MAGGING_PIPELINE_codes/9. delineage_NCLDV_bins/test/delineage_candidate.csv', 'w') as output_file:
    output_file.write('\n'.join(common_strings))
