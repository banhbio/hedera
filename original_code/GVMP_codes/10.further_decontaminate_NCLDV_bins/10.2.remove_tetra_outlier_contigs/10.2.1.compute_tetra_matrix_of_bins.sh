#require program: anvio
#input files: NCLDV bin fasta files after step 8/step 9

for FILE in /path/to/bin/*.fna;
   do
   anvi-script-compute-ani-for-fasta -f $FILE -o ${FILE%.*}_tetra_matrix.tsv -T 16 --log-file ${FILE%.*}_tetra_log --method TETRA #generate tetranucleotide matrix of NCLDV bins by anvio
   done