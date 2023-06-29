# Three Input files: 1. 'coverage_outlier_contaminated_contig_list.txt' generated from 10.2.3. ; 2.'tetra_outlier_contaminated_contig_list.txt' generated from 10.2.4. ; 3. '*new_contig_list.txt' generated from 10.1.1.

# Integrate the contaminate contigs identified from 10.1. and 10.2.
cat coverage_outlier_contaminated_contig_list.txt tetra_outlier_contaminated_contig_list.txt > contaminated_contig.list
sort contaminated_contig.list | uniq > contaminated_contig_sorted.list # the 'comm' command below needs the input file to be sorted file

# Retrieve the clean contig list in each bin (exclude the contaminate contigs)
for FILE in /path/to/bins/*new_contig_list.txt;
   do
   sort $FILE> ${FILE%.*}_sorted.txt # the 'comm' command below needs the input file to be sorted file
   comm -23 ${FILE%.*}_sorted.txt contaminated_contig_sorted.list > ${FILE%.*}_clean.list #output the lines that are unique to ${FILE%.*}_sorted.txt. 

