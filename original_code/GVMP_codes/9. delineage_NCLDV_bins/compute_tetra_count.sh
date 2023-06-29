for f in /path/to/candidate_delineage/*.fa;
   do
   id=`basename $f`
   ruby ./tetranucleotides_count_nishimura.rb $f ${id%.*}_tetra_count.tsv
   done