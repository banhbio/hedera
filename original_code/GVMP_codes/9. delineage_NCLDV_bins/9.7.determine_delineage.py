"""
If both of the two clusters satisfying all of the following features, the clusters will be the new NCLDV bins (the original NCLDV bin is determined to be delineaged): 
1) the cluster length > 40 kb (minimum GV genome size discovered in Paulo et al., 2020) ; 
2) the cluster has at least one GV single copy gene; 
3) over 50% of the contigs in the cluster have the GV characteristic score > 2; 
option criteria:
4) the species level taxonomy annotations given by CAT to the contigs are the same (species name) in the cluster but different from that in the other cluster of the bin.
"""

# Input files: 1.'*.depth.txt' generated from 9.4. ;2. 'merged_data_scored.csv' from 6.3 ; 3. '8_genes_hits_modified.csv' from 6.1

# Retrieve contig length by merge cluster contig list with the correponding '*.depth.txt', make sure 'depth.txt' and 'clusterA/B.list' are in the same directory
import os
import pandas as pd
import re

input_dir = '/Users/yue/Documents/masters_project/P002_MAGGING_PIPELINE/MAGGING_PIPELINE_codes/9. delineage_NCLDV_bins/test'
#output_dir = '/path/to/output/directory'

file_list = os.listdir(input_dir)
# Assess whether cluster A has met delineage criteria
for file_name in file_list:
    if file_name.endswith('depth.txt'):
        file_path = os.path.join(input_dir, file_name)
        
        # Modify the file_name to get the corresponding file2 path
        clusterA_name = file_name.replace('depth.txt', 'newick_clusterA.list')
        clusterA_path = os.path.join(input_dir, clusterA_name)
        df1 = pd.read_csv(clusterA_path, delimiter='\t', header=None)
        
        df2 = pd.read_csv(file_path, delimiter='\s+', header=None)
        df2 = df2.iloc[:, :2]

        # Merge df2 and df1 by their first column 
        merged_df = df2.merge(df1, left_on=0, right_on=0)

        df3= pd.read_csv("/Users/yue/Desktop/test/8_genes_hits_modified.csv", delimiter=',', header=None)
        merged_df2 = merged_df.merge(df3, left_on=0, right_on=0)
        df4= pd.read_csv("/Users/yue/Desktop/test/merged_data_scored.csv", delimiter=',', header=None)
        merged_df3 = merged_df2.merge(df4, left_on=0, right_on=0)
        merged_df3.drop_duplicates(subset=[0], inplace=True)
        
        # Rename column names
        merged_df3.columns = ['contig_id', 'contig_length', 'hallmark_gene_hits', 'max_score_group', 'CAT_data', 'score', 'query_name', 'GV_score']

       # print(merged_df3)
        # Calculate the sum of the second column
        sum_second_column = merged_df3['contig_length'].sum() > 40000

        # Check if the third column contains the specified string
        contains_string = any(merged_df3['hallmark_gene_hits'].str.contains('MCP|polB|TFIIS|VLTF3|pATPase', na=False))
        # Convert 'GV_score' column to numeric type
        merged_df3['GV_score'] = pd.to_numeric(merged_df3['GV_score'], errors='coerce')
        # Check if over 50% of the values in the last column are >= 2
        over_50_percent = (merged_df3['GV_score'] > 2).mean() > 0.5
        # Extract the string after 'family' in 'CAT_data' column if it contains 'Nucleocytoviricota'
        merged_df3['CAT_data'] = merged_df3['CAT_data'].astype(str)

        family_rows = merged_df3['CAT_data'].str.contains('family', na=False) & merged_df3['CAT_data'].str.contains('Nucleocytoviricota', na=False)
        family_df = merged_df3.loc[family_rows, 'CAT_data']

        if not family_df.empty:
           pattern = r'family\s*(\S+)'
           extracted_values = family_df.str.extract(pattern, flags=re.IGNORECASE)
           extracted_values = extracted_values.dropna()
        # Output filename with prefix matching the input filename
        output_filename = file_name.replace('depth.txt', 'clusterA_delineage_assessment.txt')
        output_path = os.path.join(input_dir, output_filename)

        # Open the output file for writing
        with open(output_path, 'w') as f:
            f.write(f"clusterB length > 40k: {sum_second_column}\n")
            f.write(f"SCG present? {contains_string}\n")
            f.write(f"Are over 50% of the contigs with GV score > 2? {over_50_percent}\n")
            if not extracted_values.empty:
                for value in extracted_values[0]:
                    f.write(f"NCLDV taxonomy annotated by CAT {value}\n")
            else:
                 f.write("no species-level annotation given by CAT")
 

# Assess whether cluster B has met delineage criteria
for file_name in file_list:
    if file_name.endswith('depth.txt'):
        file_path = os.path.join(input_dir, file_name)
        
        # Modify the file_name to get the corresponding file2 path
        clusterB_name = file_name.replace('depth.txt', 'newick_clusterB.list')
        clusterB_path = os.path.join(input_dir, clusterB_name)
        df1 = pd.read_csv(clusterB_path, delimiter='\t', header=None)
        
        df2 = pd.read_csv(file_path, delimiter='\s+', header=None)
        df2 = df2.iloc[:, :2]

        # Merge df2 and df1 by their first column 
        merged_df = df2.merge(df1, left_on=0, right_on=0)
        
        df3= pd.read_csv("/Users/yue/Desktop/test/8_genes_hits_modified.csv", delimiter=',', header=None)
        merged_df2 = merged_df.merge(df3, left_on=0, right_on=0)
    
        df4= pd.read_csv("/Users/yue/Desktop/test/merged_data_scored.csv", delimiter=',', header=None)
        merged_df3 = merged_df2.merge(df4, left_on=0, right_on=0)
        merged_df3.drop_duplicates(subset=[0], inplace=True)
        
        
        # Rename column names
        merged_df3.columns = ['contig_id', 'contig_length', 'hallmark_gene_hits', 'max_score_group', 'CAT_data', 'score', 'query_name', 'GV_score']
        
        
        # Calculate the sum of the second column
        sum_second_column = merged_df3['contig_length'].sum() > 40000

        # Check if the third column contains the specified string
        contains_string = any(merged_df3['hallmark_gene_hits'].str.contains('MCP|polB|TFIIS|VLTF3|pATPase', na=False))
        # Convert 'GV_score' column to numeric type
        merged_df3['GV_score'] = pd.to_numeric(merged_df3['GV_score'], errors='coerce')
        # Check if over 50% of the values in the last column are >= 2
        over_50_percent = (merged_df3['GV_score'] > 2).mean() > 0.5
        # Extract the string after 'family' in 'CAT_data' column if it contains 'Nucleocytoviricota'
        merged_df3['CAT_data'] = merged_df3['CAT_data'].astype(str)

        family_rows = merged_df3['CAT_data'].str.contains('family', na=False) & merged_df3['CAT_data'].str.contains('Nucleocytoviricota', na=False)
        family_df = merged_df3.loc[family_rows, 'CAT_data']

        if not family_df.empty:
           pattern = r'family\s*(\S+)'
           extracted_values = family_df.str.extract(pattern, flags=re.IGNORECASE)
           extracted_values = extracted_values.dropna()

        # Output filename with prefix matching the input filename
        output_filename = file_name.replace('depth.txt', 'clusterB_delineage_assessment.txt')
        output_path = os.path.join(input_dir, output_filename)

        # Open the output file for writing
        with open(output_path, 'w') as f:
            f.write(f"clusterB length > 40k: {sum_second_column}\n")
            f.write(f"SCG present? {contains_string}\n")
            f.write(f"Are over 50% of the contigs with GV score > 2? {over_50_percent}\n")
            if not extracted_values.empty:
                for value in extracted_values[0]:
                    f.write(f"NCLDV taxonomy annotated by CAT {value}\n")
            else:
                 f.write("no species-level annotation given by CAT")
