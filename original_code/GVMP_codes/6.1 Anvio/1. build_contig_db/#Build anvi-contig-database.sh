#Build anvi-contig-database
#Anvio Fasta file format requirement: contig title requirement-no space/'|'' etc.
#Note that this cmd do not overwrite on existing db, therefore avoid same contig db name in the same directory
module load anvio / conda activate anvio
anvi-gen-contigs-database -f contigs.fa -o contigs.db -n 'An example contigs database'

