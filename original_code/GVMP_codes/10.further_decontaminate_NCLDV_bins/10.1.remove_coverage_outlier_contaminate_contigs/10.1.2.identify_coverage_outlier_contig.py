#input files : *contig_depth.csv generated from 10.1.1.

import pandas as pd
import glob

def identify_outliers(csv_file, output_file):
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_file, header=None, names=['ID', 'Value'])

    # Calculate the first quartile (Q1), third quartile (Q3), and IQR
    Q1 = df['Value'].quantile(0.25)
    Q3 = df['Value'].quantile(0.75)
    IQR = Q3 - Q1

    # Identify outliers based on the IQR
    outliers = df[(df['Value'] < Q1 - 1.5 * IQR) | (df['Value'] > Q3 + 1.5 * IQR)]

    # Save the IDs of the outliers to the output file
    outliers['ID'].to_csv(output_file, index=False)

    # Return whether outliers were found or not
    return len(outliers) > 0

# Directory containing the files
directory = '/path/to/*contig_depth.csv/'

# Find all files ending with 'contig_depth.csv' in the directory
files = glob.glob(directory + '*contig_depth.csv')

# Initialize a flag to check if any outliers were found
outliers_found = False

# Process each file
for file in files:
    # Create the output file path by appending '_outliers.txt' to the original file name
    output_file = file.replace('.csv', '_outliers.txt')

    # Identify outliers and check if any were found
    if identify_outliers(file, output_file):
        outliers_found = True
        print(f"Outliers found in {file}")
    else:
        print(f"No outliers found in {file}")

# Print overall summary
if outliers_found:
    print("Outliers found in at least one file.")
else:
    print("No outliers found in any file.")

