# hedera
Author: Hiroki Ban (ban@kuicr.kyoto-u.ac.jp)

This is a WIP pipeline for binning NCLDV MAGs.
Original source code was provided by Ivy (Yue Fang).

## Installation
1. Install Nextflow (https://www.nextflow.io/) in your PATH.
2. prepare CAT database in data/CAT.
3. modify virarecall script to run everywhere.

## progress

1. [x] reads quality filter
2. [x] de novo assembly (contigs)
3. [x] mapping
4. [x] binning
5. [x] detect putative NCLDV bins
6. [ ] Assessment of GV characteristics of the contigs
7. [ ] discard nonNCLDV bins
8. [ ] premove cellular contigs
9. [ ] delineage NCLDV bins
10. [ ] further decontaminate NCLDV bins


## Concern
I donot care any licences yet (hmms and modefied viralrecall scripts etc.).
At least, hmm file should be replaced by us.
Please make sure before this page become public.