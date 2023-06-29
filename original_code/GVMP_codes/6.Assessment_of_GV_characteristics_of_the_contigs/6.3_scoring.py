import pandas as pd
# Specify the directory where your files are located
file_path = "/Users/yue/Desktop/test/"

# Read the merged data from the previous step
merged_df = pd.read_csv(file_path + 'merged_data.csv')

# Initialize the 'GV_score' column with 0
merged_df['GV_score'] = 0

# Assign scores based on criteria
merged_df.loc[merged_df['max_score_group'].str.contains('NCLDV', na=False), 'GV_score'] += 1
merged_df.loc[merged_df['CAT_data'].str.contains('Nucleocytoviricota', na=False), 'GV_score'] += 1
merged_df['score'] = pd.to_numeric(merged_df['score'], errors='coerce')  # Convert 'score' column to numeric
merged_df.loc[merged_df['score'] > 0, 'GV_score'] += 1
merged_df.loc[merged_df['query_name'] != 'nah', 'GV_score'] += 1

# Cap the GV_score at 4
merged_df['GV_score'] = merged_df['GV_score'].clip(upper=4)

# Save the updated DataFrame back to disk as a .csv file
merged_df.to_csv(file_path + 'merged_data_scored.csv', index=False)

# Output the cellular contig id
# Extract the contig id where 'GV_score' is '0'
cellular_contigs = merged_df.loc[merged_df['GV_score'] == 0, 'contig'].tolist()

# Write the extracted values to a text file
with open('/Users/yue/Desktop/test/cellular_contig_id.txt', 'w') as f:
    for contig in cellular_contigs:
        f.write(contig + '\n')



