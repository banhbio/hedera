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

## Pipeline Description

This pipeline is divided into several processes, recovering NCLDV MAGs from metagenomes. It primarily consists of two stages: first, the assembly and NCLDV bin screening stage, and second, the bin quality control stage.

### Stage One: Assembly and Screening

- **Quality Check of Reads:** The quality of reads is checked, with parameters customizable via options. This step can be skipped with the `--after_qc`.
- **Assembly:** Assembly is done using Megahit, which runs in meta-large mode by default. This step can be skipped with the `--from_contig` options (`--input_contigs` are required).
- **Binning:** Binning is performed using Metabat2. During this step, contigs of more than 30Kb are also recovered from the leftovers that were not binned.
- **Filtering of NCLDV Bins:** Bins (single contig) with a "Core Gene Index" greater than 5.75 are screened as NCLDVs.

### Stage Two: Quality Control of NCLDV Bins[^1]

- **Evaluation of NCLDV Features:** The obtained bins are evaluated for their likelihood of being an NCLDV bin. Each contig in the bin is scored (NCLDV score: 0-4) based on Viralrecall, VirSorter2, Contig Annotation Tool (CAT), and hmms of 149 NCLDV-characteristic genes.
- **Removal of Non-NCLDV Bins:** Bins lacking five NCLDV hallmark genes (MCP_NCLDVs, DNApolB, TFIIS, VLTF3, pATPase_all) and where over 90% of the contigs have an NCLDV score of 0 are removed from the dataset.
- **Removal of Non-NCLDV Contigs:** Contigs within bins that have an NCLDV score of 0 are removed.
- **Splitting of Chimeric Bins:** Bins that have become chimeras of multiple NCLDVs are separated. The process for this is as follows:
   1. **Identification of Candidates:** Bins with heterogeneous depth (coefficient of variation: (standard deviation / mean) of "total average depth" > 0.01) and multiple single-copy NCLDV core genes (TFIIS, VLTF3, pATPase_all, DNApolB) are considered as candidates.
   2. **Contig Clustering:** Contigs from candidate bins are clustered based on coverage and Tetranucleotide Frequency (TF).
   3. **Decision to Split:** The decision to split or not is made based on each cluster's size (40Kb), the presence or absence of single-copy NCLDV core genes, and the NCLDV score. Bins that are not divided are marked as mscg (multiple single copy genes).
- **Second Decontamination:** Outlying contigs based on coverage, TF, and NCLDV are excluded from the bin. Contigs that satisfy the following three conditions are removed from the bin:
   1. **Average Depth Outlier:** "total average depth" < Q1-1.5 \* IQR or "total average depth" > Q3+1.5 \* IQR
   2. **Principle Component Outlier:**: After applying Principle Component Analysis to the TF, if PC1 < -2.5 std. of PCA or PC1 > 2.5 std. of PCA.
   3. **NCLDV socre:** "NCLDV score" == 1

[^1]: In this context, "depth" refers to what is used for the input of Metabat2 (i.e., depth.txt).

## Concern
I donot care any licences yet (hmms and modefied viralrecall scripts etc.).
~~At least, hmm file should be replaced by us.~~
The original sources of the hmm files were confirmed.

Please make sure before this page become public.

This pipeline makes a lot of intermediate files.
CAT is so slow.
In CAT, pithovirus is not included in Nucleocytoviricota.