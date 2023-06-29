#!/bin/bash
# Retrieve candidate bins for delineage 
# Input file: 'delineage_candiate.csv' generated from 9.1.3. ; 'depth.txt' depth table of all contigs generated from 'jgi_summarize_bam_contig_depths' program
mkdir delineage_candidate_bins

input_file="delineage_candiate.csv" # PLS CHECK whether the bin id match with the format of bin id in their fasta filename

# Source directory containing the candidate bin fasta files
source_directory="/path/to/source_directory"

# Destination directory to copy the files
destination_directory="/path/to/destination_directory/delineage_candidate_bins"

# Loop over each BIN ID in the input file
while IFS= read -r id; do
    # Search for files starting with the ID and ending with '.fa' in the source directory
    files=$(find "$source_directory" -type f -name "$id*.fa")

    # Copy the matching files to the destination directory
    cp $files "$destination_directory"
done < "$input_file"

mkdir coverage_matrix # Make a directory to store the coverage matrix

for f in delineage_candidate_bins/*.fa;
   do
   # Retrieve contig list in each bin
   grep 'k141' $f > ${f%.*}.contig.list
   sed -i s'/>//'g ${f%.*}.contig.list
   awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($1 in a)print $0"\t"a[$1]}' ${f%.*}.contig.list /path/to/depth.txt > coverage_matrix/${f%.*}.depth.txt
   done
