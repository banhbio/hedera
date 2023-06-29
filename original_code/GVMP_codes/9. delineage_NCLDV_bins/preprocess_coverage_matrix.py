import os

# Input directory containing the '*depth.txt'
input_directory = "your_input_directory"

# Output directory to store the modified files
output_directory = "your_output_directory"

# Loop over each file in the input directory
for filename in os.listdir(input_directory):
    if filename.endswith("depth.txt"):
        # Construct the input file path
        input_file = os.path.join(input_directory, filename)
        
        # Construct the output file path
        output_file = os.path.join(output_directory, f"{os.path.splitext(filename)[0]}_modified.txt")
        
        # Open the input and output files
        with open(input_file, "r") as file_in, open(output_file, "w") as file_out:
            # Process each line in the input file
            for line in file_in:
                # Split the line into fields
                fields = line.strip().split("\t")
                
                # Remove the second(contig length), third(total depth), and last field (contig id)
                modified_fields = [fields[i] for i in range(len(fields)) if i not in (1, 2, len(fields)-1)]
                
                # Write the modified fields to the output file
                file_out.write("\t".join(modified_fields) + "\n")
