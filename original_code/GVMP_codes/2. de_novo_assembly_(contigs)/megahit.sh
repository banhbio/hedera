#megahit
module load megahit
#co-assembly cmd line
#integrate reads file list
ls *_1.paired.fastq.gz > 1_paired_files.list
tr '[:space:]' ',' < 1_paired_files.list > 1_paired_files.list1  #same for 2_paired_files
#coassemble_megahit.sh  submit to supercomputer machine
""
#/bin/sh

#PBS -q SDF
#PBS -N co-assemble20211231
#PBS -l select=1:ncpus=100:mem=2000gb

source /etc/profile.d/modules.sh
module load megahit/1.2.9

cd /aptmp/ideas2/../(paired_end_files_directory)
megahit (-m 1) -t 100 (--k-list 29,39,49,59,69,79,89,99,109,119,129,141) -1  (1_paired_files) -2 (2_paired_files) -o output_path 

""
qsub coassemble_megahit.sh

#individual assembly
#batch job via qsubarray
find /aptmp/../trim_reads/ERR*_1.paired.fastq.gz | parallel --dry-run "megahit (-m 0.5) -t 10 (--presets meta-sensitive) -1 {} -2 {}_2 -o /aptmp/ideas2/ivyyue/Ura_1718/coassembly/{/.}_o" > megahit_individual_assemble.qsubarray
sed -i s'/_1.paired.fastq.gz_2/_2.paired.fastq.gz/g' megahit_individual_assemble.qsubarray
sed -i s'/_1.paired.fastq_o//g' megahit_individual_assemble.qsubarray
qsubarraypbs -q SDF -l select=1:ncpus=40:mem=240gb megahit_individual_assemble.qsubarray

#filter out contigs shorter than 2500nt
module load seqkit
seqkit seq -m 2500 1718_individual_contigs.fa > 1718_individual_contigs_2500.fa