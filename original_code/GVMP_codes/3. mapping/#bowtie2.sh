#bowtie2
#index contig set
Module load bowtie2
nohup bowtie2-build 1718_individual_contigs_2500.fa 1718_individual_contigs_2500 & ; #or put on job to 'screen'
#mapping:generate SAM files
#batch jobs via qsubarray
find /aptmp/../UU*_1.paired.fastq.gz | parallel --dry-run "bowtie2 -p 5 -x /aptmp/../1718_individual_contigs_renamed_2500 -1 {} -2 {}_2 -S /aptmp/../{/.}.sam" > bowtie2.qsubarray
sed -i s'/_1.paired.fastq.gz_2/_2.paired.fastq.gz/g' bowtie2.qsubarray
sed -i s'/_1.paired.fastq.sam/.sam/g' bowtie2.qsubarray
qsubarraypbs -q APC -l select=1:ncpus=5:mem=50gb bowtie2.qsubarray

