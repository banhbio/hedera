
# Input file: NCLDV bin fna files & 'cellular_contig_id.txt'
# extract contig id list in each NCLDV bin first
for FILE in /path/to/bins/*.fna;
   do
   grep 'k141' $FILE > ${FILE%.*}_contig_id.txt #extract contig id in every NCLDV bins
   sed -i s'/>//'g ${FILE%.*}_contig_id.txt #retain contig id only
   sort ${FILE%.*}_contig_id.txt > ${FILE%.*}_contig_id_sorted.txt
   done

# retreive non-cellular contig id in each NCLDV bins
for f in *contig_id_sorted.txt;
      do
      comm -23 $f cellular_contig_id.txt > ${f%.*}_NONcellular.list #output the lines that are unique to $f. 
      done