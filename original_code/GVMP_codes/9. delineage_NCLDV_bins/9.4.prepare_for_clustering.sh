# Retrieve delineage candidate bins and compute coverage matrix
module load Python/3.8.7
bash get_contig_depths_each_bin.sh
python preprocess_depth_matrix.py 

# Compute tetranucleotide matrix
bash compute_tetra_count.sh
python normalize_tetra_count.py #make tetra matrix tsv to csv (modify this in the py script)
# Modify bin filename if you need  * Bin filename format: bin_id.additioanl.ext
#rename bin_ bin. *tetra_count_normalized.csv # make sure *depth_modified.txt and *tetra_count_normalized.csv have the same prefix e.g. 'metabat2bin.185'

mkdir newick # Store the second clustering newick results

