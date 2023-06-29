import pandas as pd

# Specify the directory where your files are located
file_path = "/Users/yue/Desktop/test/"

# Load the csv data into pandas DataFrames
viralrecall_annotation_df = pd.read_csv(file_path + 'viralrecall_annotation_fixed.csv')
final_viral_score_df = pd.read_csv(file_path + 'final-viral-score.csv')
cat_annotation_df = pd.read_csv(file_path + 'CAT_annotation_modified.csv')
hmm_hits_df = pd.read_csv(file_path + '149hmms_modified.csv')


# Fill NaN values in the first row with a NA string
final_viral_score_df.iloc[0] = final_viral_score_df.iloc[0].fillna('NA')
cat_annotation_df.iloc[0] = cat_annotation_df.iloc[0].fillna('NA')
viralrecall_annotation_df.iloc[0] = viralrecall_annotation_df.iloc[0].fillna('NA')
hmm_hits_df.iloc[0] = hmm_hits_df.iloc[0].fillna('NA')

# # Identify columns with 'k141' in their first entry
k141_col_final_viral_score = final_viral_score_df.columns[(final_viral_score_df.iloc[0].astype(str).str.contains('k141'))].tolist()[0]
k141_col_cat_annotation = cat_annotation_df.columns[(cat_annotation_df.iloc[0].astype(str).str.contains('k141'))].tolist()[0]
k141_col_viralrecall_annotation = viralrecall_annotation_df.columns[(viralrecall_annotation_df.iloc[0].astype(str).str.contains('k141'))].tolist()[0]
#k141_col_hmm_hits = hmm_hits_df.columns[(hmm_hits_df.iloc[0].astype(str).str.contains('k141'))].tolist()[0]
#print(hmm_hits_df.iloc[0])
# Extract the required columns and rename them
final_viral_score_df = final_viral_score_df[[k141_col_final_viral_score, 'max_score_group']]
final_viral_score_df.columns = ['contig', 'max_score_group']
final_viral_score_df['contig'] = final_viral_score_df['contig'].str.replace('full', '')
final_viral_score_df['contig'] = final_viral_score_df['contig'].str.replace(r'\|\|', '')

cat_annotation_df = cat_annotation_df.rename(columns={k141_col_cat_annotation: 'contig'})

viralrecall_annotation_df = viralrecall_annotation_df[[k141_col_viralrecall_annotation, 'score']]
viralrecall_annotation_df.columns = ['contig', 'score']

# Process the 'eight_genes_hits_df' DataFrame
hmm_hits_df['contig'] = hmm_hits_df['target_name']
hmm_hits_df = hmm_hits_df[['contig', 'query_name']]

# Merge all the dataframes on 'k141'
merged_df = pd.merge(final_viral_score_df, cat_annotation_df, on='contig', how='outer')
merged_df = pd.merge(merged_df, viralrecall_annotation_df, on='contig', how='outer')
merged_df = pd.merge(merged_df, hmm_hits_df, on='contig', how='outer')

# Fill NA values with 'nah'
merged_df = merged_df.fillna('nah')

# Remove rows where 'k141' is not present in the 'k141'
merged_df = merged_df[merged_df['contig'].str.contains('k141')]

# Save the merged DataFrame back to disk as a .csv file
merged_df.to_csv(file_path + 'merged_data.csv', index=False)
