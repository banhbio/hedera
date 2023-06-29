# Input files: '*.out' generated from 6.0.5.
#!/bin/bash

input_directory="/path/to/*.out"
output_file="/path/to/8_genes_hits.csv"

# Create the output file with column headers
echo -e "contig_id,gene" > "$output_file"

# Loop through the input files
for file in "$input_directory"/*.out; do
    # Extract the contig id and gene columns
    awk -F ' +' '!/^#/ && match($1, /k141_[^_]+/) { seqname = substr($1, RSTART, RLENGTH); gene = $3; print seqname "," gene }' "$file" >> "$output_file"

done

echo "Output file has been created."
