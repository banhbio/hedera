#generate a depth file from sorted BAM file
module load MetaBAT2 (?)
jgi_summarize_bam_contig_depths --outputDepth depth.txt *.sorted.bam
