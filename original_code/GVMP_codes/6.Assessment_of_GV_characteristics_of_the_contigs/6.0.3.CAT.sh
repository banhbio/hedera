module load CAT
module load Python/3.8.7
module load diamond
module load prodigal
CAT contigs -c xxx.fa -d /aptmp/ideas2/ivyyue/CAT/CAT_pack/CAT_prepare_20210107/2021-01-07_CAT_database/ -t /aptmp/ideas2/ivyyue/CAT/CAT_pack/CAT_prepare_20210107/2021-01-07_taxonomy
#add names to taxonomy ids
CAT add_names -i out.CAT.contig2classification.txt -o individual_putative_NCLDV_bins_contigs_CAT_annotation.txt -t /aptmp/ideas2/ivyyue/CAT/CAT_pack/CAT_prepare_20210107/2021-01-07_taxonomy
#replace space to | in annotation file
cut -f2- ./143_contigs_CAT_annotation.txt | sed -r s'/[[:space:]]/|/g' > ./143_contigs_CAT_annotation2.txt
cut -f1 ./143_contigs_CAT_annotation.txt > ./143_contigs_CAT_annotation1.txt
paste -d, ./143_contigs_CAT_annotation1.txt ./143_contigs_CAT_annotation2.txt > ./143_contigs_CAT_annotation_modified.csv