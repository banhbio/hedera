# Input file: fna files of the putative NCLDV bins

for FILE in /path/to/putative/NCLDV/bins/*.fna;
   do
   grep 'k141' $FILE > ${FILE%.*}.tmp1 #extract contig ids in each putative NCLDV bins
   sed -i s'/>//'g ${FILE%.*}.tmp1 #retain contig id only
   awk '{print $0 "," FILENAME}' ${FILE%.*}.tmp1 > ${FILE%.*}.tmp2
   cat ${FILE%.*}.tmp2 > all_putative_bins_contig_list.csv
   done

rm -f /path/to/putative/NCLDV/bins/*.tmp1
rm -f /path/to/putative/NCLDV/bins/*.tmp2 
sed -i s'/\.tmp1//'g all_putative_bins_contig_list.csv #retain bin id only