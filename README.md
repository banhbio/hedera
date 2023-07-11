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
6. [x] assessment of GV characteristics of the contigs
7. [x] discard nonNCLDV bins
8. [x] remove cellular contigs
9. [x] delineage NCLDV bins (except CAT)
10. [x] further decontaminate NCLDV bins



## Concern
I donot care any licences yet (hmms and modefied viralrecall scripts etc.).
~~At least, hmm file should be replaced by us.~~
The original sources of the hmm files were confirmed.

Please make sure before this page become public.

This pipeline makes a lot of intermediate files.
CAT is so slow.
In CAT, pithovirus is not included in Nucleocytoviricota.