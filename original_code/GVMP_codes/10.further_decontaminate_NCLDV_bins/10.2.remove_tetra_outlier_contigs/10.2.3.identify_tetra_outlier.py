#identify tetra outlier which show PC1 value < −2.5 s.d. or > 2.5 s.d. (Nishimura et al., 2022). 
import os
import pandas as pd
import numpy as np

# Set the input and output directories
input_directory = "/path/to/input/directory"
output_directory = "/path/to/output/directory"

# Get a list of input files in the directory
input_files = os.listdir(input_directory)

# Loop through the input files
for input_file in input_files:
    # Read the input file into a DataFrame
    data = pd.read_csv(os.path.join(input_directory, input_file))
    
    # Calculate the standard deviation of the PC1_Value column
    pc1_std = np.std(data['PC1_Value'])
    
    # Find the outliers based on the specified criteria
    outliers = (data['PC1_Value'] < -2.5 * pc1_std) | (data['PC1_Value'] > 2.5 * pc1_std)
    
    # Check if any outliers were found
    if outliers.any():
        # Filter the output dataframe to include only the outliers
        outlier_df = data[outliers]
        
        # Remove the double quotation marks from the Contig_ID column
        outlier_df["Contig_ID"] = outlier_df["Contig_ID"].str.replace('"', "")
        
        # Create the output file name
        output_file = os.path.join(output_directory, os.path.splitext(input_file)[0] + "_outliers.csv")
        
        # Write the outliers to the output file
        outlier_df["Contig_ID"].to_csv(output_file, index=False)
        
        # Print a message indicating the output file name
        print("Outliers found in", input_file, "and written to", output_file)
    else:
        # No outliers found, print a message indicating the input file
        print("No outliers found in", input_file)
