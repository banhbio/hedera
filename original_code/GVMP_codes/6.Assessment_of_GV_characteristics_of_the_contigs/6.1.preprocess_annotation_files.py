import re
import pandas as pd
# Specify the directory where your files are located
file_path = "/Users/yue/Desktop/test/"

# List all the files that need to be pre-processed
files = ['final-viral-score.tsv', 'CAT_annotation.txt', 'viralrecall_annotation.summary.tsv']

# For each file
for file in files:
    # Open the file and read the data
    with open(file_path + file, 'r') as f:
        lines = f.readlines()

    # For each line, replace all sequences of spaces or tabs with a single comma
    lines = [re.sub('\s+', ',', line.strip()) for line in lines]

    # Write the processed lines to a new .csv file
    with open(file_path + file.split('.')[0] + '.csv', 'w') as f:
        f.write('\n'.join(lines))



## pre-process CAT_annotation.csv since it got inconsistency and complex text in lines
output_path = "/Users/yue/Desktop/test/CAT_annotation_modified.csv"

k141_ids = []
root_data = []

with open(file_path + 'CAT_annotation.csv', 'r') as file:
    for line in file:
        # Skip lines that don't start with 'k141'
        if not line.startswith('k141'):
            continue
        
        # Split the line into comma-separated fields
        fields = line.strip().split(',')

        # The first field is the 'k141' identifier
        k141_ids.append(fields[0])

        # Look for 'root' field, get its index and everything after that
        if 'root' in fields:
            root_index = fields.index('root')
            root_data.append(','.join(fields[root_index:]))
        else:
            root_data.append('N/A')

df = pd.DataFrame({'k141': k141_ids, 'CAT_data': root_data})

# Write the DataFrame to a CSV file
df.to_csv(output_path, index=False)

# Fix broken line issues in viralrecall_annotation.csv
with open(file_path + 'viralrecall_annotation.csv', 'r') as file:
    lines = file.readlines()

# Fix the broken lines
fixed_lines = []
current_line = ""
for line in lines:
    if "k141_" in line:
        if current_line:
            fixed_lines.append(current_line)
        current_line = line.strip()
    else:
        current_line += line.strip()
fixed_lines.append(current_line)  # don't forget the last line

# Save the fixed data to a new file
with open(file_path + 'viralrecall_annotation_fixed.csv', 'w') as file:
    file.write("\n".join(fixed_lines))

### pre-process 149hmms.out to deal with duplicate annotation
# Specify the file path
path = "/Users/yue/Desktop/test/149hmms.out"

# Read the file skipping the header rows
df = pd.read_csv(path, skiprows=3, delim_whitespace=True, usecols=[0, 2], names=['target_name', 'query_name'])

# Retain content before the second '_'
df['target_name'] = df['target_name'].str.split('_', n=2).str[:2].str.join('_')

# Group by first column and merge the second column
df = df.groupby('target_name')['query_name'].apply(':'.join).reset_index()

# Save the DataFrame to a new CSV file
output_file = "/Users/yue/Desktop/test/149hmms_modified.csv"
df.to_csv(output_file, index=False)

### pre-process 8_genes_hits.csv to deal with duplicate annotation
# Specify the file path
path = "/Users/yue/Desktop/test/8_genes_hits.csv"

# Read the file skipping the header rows
df = pd.read_csv(path)

# Retain content before the second '_'
# df['target_name'] = df['target_name'].str.split('_', n=2).str[:2].str.join('_')

# Group by first column and merge the second column
df = df.groupby('seqname')['gene'].apply(':'.join).reset_index()

# Save the DataFrame to a new CSV file
output_file = "/Users/yue/Desktop/test/8_genes_hits_modified.csv"
df.to_csv(output_file, index=False)


















