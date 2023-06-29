# PLEASE CHECK whether contig id in *tetra_matrix_PC1_results_outliers.csv match with the contig id in 'merged_data_scored.csv' generated from 6.3.script

# concatenate all '*tetra_matrix_PC1_results_outliers.csv'
cat *tetra_matrix_PC1_results_outliers.csv > tetra_outlier_concat.txt # concatenate all the contig list files into one file
awk 'BEGIN{FS=","} NR==FNR{a[$1]=$6;next}NR>FNR{if($1 in a)print $1","$6 a[$1]}' tetra_outlier_concat.txt /path/to/merged_data_scored.csv >  tetra_outlier_scores.csv

# now extract the contig id with GV score of '1' which are considered as contaminate contigs
awk -F',' '$2 == 1 {print $1}' tetra_outlier_scores.csv > tetra_outlier_contaminated_contig_list.txt