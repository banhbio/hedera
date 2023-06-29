cd /aptmp/ideas2/ivyyue/viralrecall #go to viralrecall.py directory
module load hmmer
module load prodigal
module load Biopython
module load python/3.8.7  /3.7.5
screen
python viralrecall.py -i /path/contigs.fa -p /output/path/prefix  -c(contig level) -t 32 &