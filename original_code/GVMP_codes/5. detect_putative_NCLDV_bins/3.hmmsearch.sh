mkdir all_bins_hmm_result
find /aptmp/lingjie/001_Research/06_Uranouchi/10_binning/hmm_file/hallmark_hmm/*.hmm | parallel --dry-run "hmmsearch --cpu 10 --notextw --tblout /outputpath/{/.}.out {} /*/all_bins_cds.faa" > hmmsearch.qsubarray
module load hmmer
qsubarraypbs -q APC -l select=1:ncpus=10:mem=50gb hmmsearch.qsubarray
#Output: /path/all_individual_bins_hmm_result
#modify_ids..
#sed -i s'/metabat/_metabat/g' *.out
#sed -i s'/.genes_k141_[[:digit:]]+//g' *.out
#NOTICE:bin filename should be same as the hits id in hmm_results.out (exclude suffix)