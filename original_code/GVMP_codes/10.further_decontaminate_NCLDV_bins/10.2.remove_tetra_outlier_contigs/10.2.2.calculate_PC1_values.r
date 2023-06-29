# Load required library
library(stats)
library(tools)

# Set the input and output directories
input_directory <- "/Users/yue/Desktop/test/"  
output_directory <- "/Users/yue/Desktop/test/"  

# Get a list of input files in the directory
input_files <- list.files(input_directory, pattern =".*tetra_matrix\\.tsv$", full.names = TRUE)

# Loop through the input files
for (input_file in input_files) {
  # Read the matrix from the input file
  data <- read.csv(input_file, sep = "\t", header = TRUE)
  
  # Perform PCA
  pca <- prcomp(data[,2:ncol(data)], scale=F)
  
  # Extract PC1 values for each contig
  contig_ids <- data[, 1] # store contig ids
  pc1_values <- pca$x[, 1]
  
  # Create the output file name
  output_file <- file_path_sans_ext(basename(input_file))
  output_file <- paste0(output_file, "_PC1_results.csv")
  output_file <- file.path(output_directory, output_file)
  
  # Combine contig IDs and PC1 values into a data frame
  output_df <- data.frame(Contig_ID = contig_ids, PC1_Value = pc1_values)
  
  # Output the data frame to a CSV file
  write.csv(output_df, file = output_file, row.names = FALSE)
  
  # Print a message indicating the output file name
  cat("Output file", output_file, "has been created.\n")
}







