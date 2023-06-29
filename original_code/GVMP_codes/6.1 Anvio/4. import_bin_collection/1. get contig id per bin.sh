#I.get contig id per bin
""

#!/bin/bash
# A script to convert binning results for anvio input

FILES=$(find *.fa)
for f in $FILES; do
 NAME=$(basename $f .fa)
 grep ">" $f | sed 's/>//' | sed -e "s/$/\t$NAME/" | sed 's/\./_/' >> contigs_per_bin.tsv
done

""