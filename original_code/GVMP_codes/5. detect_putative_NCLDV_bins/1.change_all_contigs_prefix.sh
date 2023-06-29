#contig id 前加上bin id
for f in ./*.fa;
    do
    id=`basename ${f%.*}`
    sed -i "s/k141/$id\_k141/" $f
    done
     #rename metabat "$sid"_metabat $dir*.fa

#去掉‘#’后的str 
sed -i s'/#.*//g' all_bins_cds.faa