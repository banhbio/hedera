# Iterately run second clustering for the delineage candidate bins based on their coverage matrix and tetra matrix
for f in coverage_matrix/*depth_modified.txt;
    do
    basename=$(basename "$f" .txt)
    basename=${basename%.*}
   # echo "$basename" # To check whether this output the correct prefix for second_clustering.py
    binid=`echo "$basename"` # make sure *depth_modified.txt and *tetra_count_normalized.csv have the same prefix e.g. 'metabat2bin.185'
    python second_clustering.py $binid
    done
    