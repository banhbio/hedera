# Three input file:  1. new version of NCLDV bin fna files generated after step 8 ; 2. *contig_list.txt generated from 8.1.get_cellular_contigs_list_each_bin.sh ;3. 'depth.txt' depth table of all contigs generated from 'jgi_summarize_bam_contig_depths' program

for FILE in /path/to/bins/*.fna;
   do
   grep 'k141' $FILE > ${FILE%.*}_new_contig_list.txt #retrieve contig id in every NCLDV bins (after 1st step decontamination and delineage)
   sed -i s'/>//'g ${FILE%.*}_new_contig_list.txt 
   awk 'NR==FNR{a[$1]=$3;next}NR>FNR{if($1 in a)print $1","$3 a[$1]}' ${FILE%.*}_new_contig_list.txt /path/to/depth.txt > ${FILE%.*}_contig_depth.csv # retrieve the total depth of each contig in the bins
   done

# sample content of output file
# k141_4948495,46.4459
# k141_2480350,63.8615
# k141_5010737,47.2984
# k141_2524980,47.4765
