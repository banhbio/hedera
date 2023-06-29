# Remove contaminate contigs from the bins using seqkit
# Require tool kit: seqkit

for bin in /path/to/bins/*.fa;
    do
    a=`basename $bin`
    t=${FILE%.*}_clean.list
    seqkit faidx $bin -l $t > ${a%.*}_final.fa #extract the clean contigs in the NCLDV bins and make FINAL decontaminated versions of the NCLDV bins
    done