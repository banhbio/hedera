# Input files: NCLDV bin fna files & files named '*contig_id_sorted_NONcellular.list'
# Require tool kit: seqkit

   for bin in /path/to/bins/*.fa;
        do
        a=`basename $bin`
        t=${a%.*}_contig_id_sorted_NONcellular.list
        seqkit faidx $bin -l $t > ${a%.*}_removed_cellular.fa #extract the non-cellular contigs in the NCLDV bins and make new versions of the NCLDV bins
        done

