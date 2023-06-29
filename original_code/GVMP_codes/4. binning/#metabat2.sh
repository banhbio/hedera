#metabat2
module load MetaBAT2
#metabat2.sh
""

#!/bin/bash

#PBS -q APC
#PBS -N metabat2
#PBS -l select=1:ncpus=20:mem=200gb
source /etc/profile.d/modules.sh
module load MetaBAT2

cd /aptmp/ideas2/ivyyue/R9_coassem_metabat2
metabat2 -i contigs.fa -a depth.txt -v -o R9/metabat2bin #output bin preffix

""
qsub metabat2.sh