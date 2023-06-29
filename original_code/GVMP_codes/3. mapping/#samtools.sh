#samtools
#convert SAM files to BAM files
module load samtools
#batch jobs via qsubarray
find /aptmp/../*.sam | parallel --dry-run "samtools view -b -@ 5 {} > {.}.bam" > sam_to_bam.qsubarray
qsubarraypbs -q SDF -l select=1:ncpus=40:mem=240gb sam_to_bam.qsubarray
#sort and index BAM files
find /aptmp/../*.bam | parallel --dry-run "samtools sort {} -o {.}.sorted.bam -@ 5 -m 50G" > sort_bam.qsubarray
qsubarraypbs -q APC -l select=1:ncpus=6:mem=50gb sort_bam.qsubarray

find /aptmp/../*.sorted.bam | parallel --dry-run "samtools index {} -@ 5" > index_bam.qsubarray
qsubarraypbs -q SDF -l select=1:ncpus=40:mem=240gb index_bam.qsubarray