# Calculate coefficient of variation (CV) of the contigs' total depth in each bin to identify the NCLDV bins with uneven coverage distribution
# CV = (Standard Deviation / Mean of contigs' total depth) * 100
# Input files: *contig_depth.csv generated from 9.1.1.
import pandas as pd
import glob

# Specify the directory where your files are located
file_directory = "/Users/yue/Documents/masters_project/P002_MAGGING_PIPELINE/MAGGING_PIPELINE_codes/10.further_decontaminate_NCLDV_bins/10.1.remove_coverage_outlier_contaminate_contigs/test/"
output_directory = "/Users/yue/Documents/masters_project/P002_MAGGING_PIPELINE/MAGGING_PIPELINE_codes/9. delineage_NCLDV_bins/test/"
# Create an empty dataframe to store the CV values and filenames
cv_dataframe = pd.DataFrame(columns=['BIN_ID', 'CV'])

# Get the list of file paths for the files in the directory
file_paths = glob.glob(file_directory + "*_contig_depth.csv")

# Iterate over each file
for file_path in file_paths:
    # Read the file into a dataframe
    df = pd.read_csv(file_path, header=None, names=['Column1', 'Column2'])
    
    # Calculate the coefficient of variation (CV) for the second column
    cv = df['Column2'].std() / df['Column2'].mean() * 100
    
    # Get the filename from the file path and remove '_contig_depth.csv' part
    bin_id = file_path.split("/")[-1].replace("_contig_depth.csv", "")
    
    # Append the CV value and filename to the dataframe
    cv_dataframe = cv_dataframe.append({'BIN_ID': bin_id, 'CV': cv}, ignore_index=True)


# Save the dataframe to a new CSV file
cv_dataframe.to_csv(output_directory + "CV_per_bin.csv", index=False)

