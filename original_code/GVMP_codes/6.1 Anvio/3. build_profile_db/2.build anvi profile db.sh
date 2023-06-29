#build anvi profile db
#note that contig id in the bam files should be same as the contig files input to contig.db
#anvi-profile is stored automatically in the bam file directory
anvi-profile -i SAMPLE-01.bam (sorted Bam file) -c contigs.db -T 8  (threads) #with BAM file index-filename end in‘.bat’ in the same directory
#merge anvi profile db if have multiple samples
anvi-merge */PROFILE.db -o SAMPLES-MERGED -c contigs.db
