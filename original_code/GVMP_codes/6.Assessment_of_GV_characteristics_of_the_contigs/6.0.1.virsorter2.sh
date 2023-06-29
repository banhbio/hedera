#!/bin/bash

#PBS -q SDF
#PBS -N virsorter2
#PBS -l select=1:ncpus=32:mem=240gb

source /etc/profile.d/modules.sh
module purge
module load VirSorter/2.2.1

cd /aptmp/ideas2/ivyyue/Ura_1718_42_samples/binning_MetaBAT2/putative_NCLDV_filter/42_samples_putative_NCLDV_bins/contig_annotations/
virsorter run -w output_file.out -i concat_contigs.fa --include-groups "dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae" -j 32 all