#input files: 1. 'merged_data_scored.csv' generated from 6.3.script ; 2. '*contig_depth_outliers.txt' generated from 10.1.2. script


cat *contig_depth_outliers.txt > coverage_outlier_concat.txt #concatenate all the coverage outlier contig id into one file
awk 'BEGIN{FS=","} NR==FNR{a[$1]=$6;next}NR>FNR{if($1 in a)print $1","$6 a[$1]}' coverage_outlier_concat.txt /path/to/merged_data_scored.csv >  coverage_outlier_scores.csv

# sample content of the output 'coverage_outlier_scores.csv':
# k141_13219451,1
# k141_10997165,1
# k141_12653549,2
# k141_13712832,1

#now extract the contig id with GV score of '1' which are considered as contaminate contigs
awk -F',' '$2 == 1 {print $1}' coverage_outlier_scores.csv > coverage_outlier_contaminated_contig_list.txt