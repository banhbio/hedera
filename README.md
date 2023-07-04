# hedera
Author: Hiroki Ban (ban@kuicr.kyoto-u.ac.jp)

This is a WIP pipeline for binning NCLDV MAGs.
Original source code was provided by Ivy (Yue Fang).

## Installation
1. Install Nextflow (https://www.nextflow.io/) in your PATH.
2. prepare CAT database in data/CAT.
3. modify virarecall script to run everywhere.

## progress

1. [x] reads_quality_filter
2. [x] de_novo_assembly_(contigs)
3. [x] mapping
4. [x] binning
5. [x] detect_putative_NCLDV_bins
6. [ ] Assessment_of_GV_characteristics_of_the_contigs
7. [ ] discard_nonNCLDV_bins
8. [ ] premove_cellular_contigs
9. [ ] delineage_NCLDV_bins
10. [ ] further_decontaminate_NCLDV_bins


## Concern
I donot care any licences yet (hmms and modefied viralrecall scripts etc.).
At least, hmm file should be replaced by us.
Please make sure before this page become public.