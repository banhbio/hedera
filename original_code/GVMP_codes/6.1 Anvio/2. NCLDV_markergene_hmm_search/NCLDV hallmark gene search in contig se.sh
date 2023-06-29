#NCLDV hallmark gene search in contig set
#First prepare NCLDV marker gene hmm files
/aptmp/ideas2/ivyyue/anvio/anvio_bams/NCLDV_coregene_hmm_tom/
/aptmp/ideas2/ivyyue/anvio/anvio_bams/NCLDV_coregene_hmm_tom/NCLDV_coregene_hmm_tom_faylward_bitscore/ #both in supercomputer
#cmd line
anvi-run-hmms -c CONTIGS.db -H HMM_DNApolB/ -T 4 #-H DIY hmm profile source; input the whole directory
