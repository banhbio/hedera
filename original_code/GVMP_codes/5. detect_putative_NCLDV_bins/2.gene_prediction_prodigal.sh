#run prodigal for each bin
module load prodigal
mkdir ./genes_prodigal/
find ./*.fa | parallel --dry-run "prodigal -p single -a {//}/genes_prodigal/{/.}.genes.faa -d {//}/genes_prodigal/{/.}.genes.fna -f gff -o {//}/genes_prodigal/{/.}.genes.gff -i {}" > ./prodigal.qsubarray
#after prodigal finish
cat ./*.faa > all_bins_cds.faa
#去掉‘#’后的str 
sed -i s'/#.*//g' all_bins_cds.faa