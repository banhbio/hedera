#!/usr/bin/env nextflow

/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl=2

version = 'wip-0.1'

params.help = false
params.resume = false
params.cpu = 1
params.out="$baseDir/out"

if (params.help) {
    log.info """
    hedera: This is a pipleline to assemble and binning NCLDV MAGs from metagenomic data.
    version ${version}
    """
}


/* coregene settings */
params.conserved_20_NCVOG_hmm="$baseDir/data/hmm/NCVOG/conserved_20_NCVOG.hmm"
params.conserved_20_NCVOGs="NCVOG0022,NCVOG0023,NCVOG0037,NCVOG0038,NCVOG0052,NCVOG0076,NCVOG0236,NCVOG0249,NCVOG0261,NCVOG0262,NCVOG0271,NCVOG0272,NCVOG0273,NCVOG0274,NCVOG0276,NCVOG1060,NCVOG1117,NCVOG1127,NCVOG1164,NCVOG1353"
params.conserved_20_NCVOG_weights="0.9,1.1,0.5,1.1,0.9,1,0.8,1,0.7,1,0.9,1,0.8,0.9,0.8,0.6,0.7,0.4,1,0.8"
params.core_gene_index=5.75

/* QC settings*/
params.virsorter_groups="dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae"
params.CAT_DB="$baseDir/data/CAT/CAT_prepare_20210107/2021-01-07_CAT_database"
params.CAT_Taxonomy="$baseDir/data/CAT/CAT_prepare_20210107/2021-01-07_taxonomy"
params.additional_149_NCVOG_hmm="$baseDir/data/hmm/NCLDV_VIRUS_149_hmm/genes.hmm"
/*
params.hallmark_hmm="$baseDir/data/hmm/hallmark/hallmark.hmm"
params.hallmark_genes="DNApolB,MCP_all,pATPase_all,Primase_all,RNAP-a_all,RNAP-b_all,TFIIS,VLTF3"
params.hallmark_hmm_score_threshold="150,80,80,80,200,200,100,80"
*/
params.hallmark_hmm="$baseDir/data/hmm/hallmark/tmp/hallmark_without_MCP.hmm"
params.hallmark_genes="DNApolB,pATPase_all,Primase_all,RNAP-a_all,RNAP-b_all,TFIIS,VLTF3"
params.hallmark_hmm_score_threshold="150,80,80,200,200,100,80"

log.info"""
"""


workflow {
    /*01 assembly and binning*/
                                
    read_ch = Channel.fromFilePairs("${params.rawread_dir}/*_{1,2}.fastq.gz", flat:true)

    fastp(
        read_ch
    )

    megahit(
        fastp.out.read
    )
    
    forward_reads_list = fastp.out.read.map{it[1]}.toList()
    backward_reads_list = fastp.out.read.map{it[2]}.toList()

    coverm(
        megahit.out.contig,
        forward_reads_list,
        backward_reads_list
    )

    metabat2_input_ch =  megahit.out.contig.combine(coverm.out.bam, by: 0)

    metabat2(
        metabat2_input_ch
    )

    /*01 finnished*/

    /*02 detect putative NCLDV bins */

    /* collect all bins */
    bin_ch = metabat2.out.bins.flatten()

    rename_bin_header(
        bin_ch
    )

    prodigal(
        rename_bin_header.out.bin
    )

    hmmsearch_with_NCVOGs(
        prodigal.out.faa
    )

    classifier_input_ch = rename_bin_header.out.bin.combine(hmmsearch_with_NCVOGs.out.tblout, by: 0)
    
    classify_NCLDV_bin(
        classifier_input_ch
    )

    classify_NCLDV_bin.out.bin /* filter putative NCLDV bins */
                .branch {
                    ncldv: it[2].readLines().last().split('\t').last().toFloat() > params.core_gene_index
                    all: true
                }
                .set {result}

    putative_ncldv_bin_and_prot_ch = result.ncldv.map{[it[0], it[1]]}.combine(prodigal.out.faa, by: 0)
    classifier_result_list = result.all.map{it[2]}.toList()
    
    summarize_NCVOG_results(classifier_result_list)
    /*02 finnished */

    /*03 aasess putative NCLDV bins*/
    putative_ncldv_bin_ch = putative_ncldv_bin_and_prot_ch.map{[it[0], it[1]]}
    putative_ncldv_prot_ch = putative_ncldv_bin_and_prot_ch.map{[it[0], it[2]]}

    viralrecall(
        putative_ncldv_bin_ch
    )

    virsorter2(
        putative_ncldv_bin_ch
    )

    CAT(
        putative_ncldv_bin_ch
    )

    hmmsearch_with_NCLDV_VIRUS_149_hmm(
        putative_ncldv_prot_ch
    )

    assessment_results_ch = viralrecall.out.tsv.combine(virsorter2.out.tsv, by: 0)
                                               .combine(CAT.out.txt, by: 0)
                                               .combine(hmmsearch_with_NCLDV_VIRUS_149_hmm.out.tblout, by: 0)
/*    
    summarize_assesssment(
        assessment_results_ch
    )
*/
    hmmsearch_with_hallmark_genes(
        putative_ncldv_prot_ch
    )

    summarize_detected_hallmark_genes(
        hmmsearch_with_hallmark_genes.out.tblout
    )
    /*03 finnished */
}

/*
WIP: fastp
WIP: -3 -W 6 -M 30 -q 20 -u 50 -n 0 -p -l 50
*/
process fastp {
    publishDir "${params.out}/quality_check", mode: 'symlink'

    input:
    tuple val(id), path("seq1.fq.gz"), path("seq2.fq.gz")

    output:
    tuple val("${id}"), path("${id}_forward_qcd.fq.gz"), path("${id}_backward_qcd.fq.gz"), emit: 'read'
    path("${id}_fastp_report.html")
    path("${id}_fastp_report.json")

    script:
    """
    fastp -i seq1.fq.gz -I seq2.fq.gz -o ${id}_forward_qcd.fq.gz -O ${id}_backward_qcd.fq.gz -w ${task.cpus} -h ${id}_fastp_report.html -j ${id}_fastp_report.json -3 -W 6 -M 30 -q 20 -u 50 -n 0 -p -l 50
    """
}

/*
what megahit kmer is usualy used?
WIP: --k-min 21 --k-max 141

Filter out contigs shorter than 2500nt
*/
process megahit {
    publishDir "${params.out}/assembly", mode: 'symlink'

    input:
    tuple val(id), path(forward), path(backward)

    output:
    tuple val(id), path("megahit/${id}.megahit.contigs.fa"), emit: 'contig'

    script:
    mem = (task.memory =~ /(\d+).GB/)[0][1] 
    """
    megahit -1 ${forward} -2 ${backward} -m 0.3 --k-min 21 --k-max 141 --k-step 12 -t ${task.cpus} -o megahit --out-prefix ${id}.megahit --min-contig-len 2500
    """
}

/*
replaced bowtie2 and samtools with coverM, because it's simpler
*/
process coverm {
    publishDir "${params.out}/coverm", mode: 'symlink'
    input:
    tuple val(id), path(contig)
    path("seq1/*")
    path("seq2/*")

    output:
    tuple val(id), path("${id}/*"), emit: 'bam'

    script:
    """
    coverm make -r ${contig} -1 seq1/* -2 seq2/* -o ${id} -t ${task.cpus}
    """
}

process metabat2 {
    publishDir "${params.out}/metabat2", mode: 'symlink'

    input:
    tuple val(id), path(contig), path("bam/*")

    output:
    path("${id}.metabat2bin/*"), emit: 'bins'
    path("${id}.depth.txt"), emit: 'depth'

    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth ${id}.depth.txt bam/*.bam
    metabat2 -i ${contig} -a ${id}.depth.txt -v -o ${id}.metabat2bin/${id}.matabat2bin
    """
}

process rename_bin_header {
    publishDir "${params.out}/bins/fasta", mode: 'symlink'

    input:
    path(bin)

    output:
    tuple val(id), path("${id}.fasta"), emit: 'bin'

    script:
    /* define new header */
    id=bin.getBaseName().replaceAll(/\./,"_")
    """
    seqkit replace -p "k141" -r '${id}_k141' ${bin} > ${id}.fasta 
    """
}

/* replace with prodigal-gv? */
process prodigal {
    publishDir "${params.out}/bins/prodigal", mode: 'symlink'

    input:
    tuple val(id), path(bin)

    output:
    tuple val(id), path("${id}.genes.faa"), emit: 'faa'

    script:
    """
    prodigal -i ${bin} -p single -a ${id}.genes.faa -d ${id}.genes.fna -f gff -o ${id}.genes.gff
    """
}

process hmmsearch_with_NCVOGs {
    publishDir "${params.out}/NCVOG/hmm", mode: 'symlink'

    input:
    tuple val(id), path(faa)

    output:
    tuple val(id), path("${id}.NCVOG.tblout"), emit: 'tblout'

    script:
    """
    hmmsearch --cpu ${task.cpus} -E 1e-03 --notextw --tblout ${id}.NCVOG.tblout ${params.conserved_20_NCVOG_hmm} ${faa}
    """
}

process classify_NCLDV_bin {
    input:
    tuple val(id), path(bin), path(tblout)

    output:
    tuple val(id), path(bin), path("${id}.20NCVOG_weight.tsv"), emit: 'bin'

    script:
    """
    python ${baseDir}/bin/classify_bin_20NCVOG.py -f ${bin} -t ${tblout} -n ${params.conserved_20_NCVOGs} -w ${params.conserved_20_NCVOG_weights} -o ${id}.20NCVOG_weight.tsv 
    """
}

process summarize_NCVOG_results {
    publishDir "${params.out}/NCVOG", mode: 'symlink'

    input:
    path("table/*")
    
    output:
    path("20NCVOG_weight.tsv")

    script:
    """
    python ${baseDir}/bin/summarize_table.py -i ./table -p .20NCVOG_weight.tsv -o 20NCVOG_weight.tsv
    """
}

process viralrecall {
    publishDir "${params.out}/assessment/viralrecall", mode: 'symlink'

    input:
    tuple val(id), path(bin)
    
    output:
    tuple val(id), path("${id}.summary.tsv"), emit:'tsv'

    script:
    """
    viralrecall.py -i ${bin} -p ${id} -c -t ${task.cpus}
    mv ${id}/${id}.summary.tsv .
    """
}

process virsorter2 {
    publishDir "${params.out}/assessment/virsorter2", mode: 'symlink'

    input:
    tuple val(id), path(bin)
    
    output:
    tuple val(id), path("${id}.virsorter_score.tsv"), emit:'tsv'

    script:
    """
    virsorter run -w ${id} -i ${bin} --include-groups ${params.virsorter_groups} -j ${task.cpus}
    mv ${id}/final-viral-score.tsv ${id}.virsorter_score.tsv
    """
}

process CAT {
    publishDir "${params.out}/assessment/CAT", mode: 'symlink'

    input:
    tuple val(id), path(bin)
    
    output:
    tuple val(id), path("${id}.CAT_annotation.txt"), emit:'txt'

    script:
    """
    CAT contigs -c ${bin} -o ${id} -d ${params.CAT_DB} -t ${params.CAT_Taxonomy} -n ${task.cpus}
    CAT add_names -i ${id}.contig2classification.txt -o ${id}.CAT_annotation.txt -t ${params.CAT_Taxonomy}
    """
}

/* Is it correct threshold? */
/* the ivy's docs said 1e-10 */
/* original code said 1e-50 */
process hmmsearch_with_NCLDV_VIRUS_149_hmm {
    publishDir "${params.out}/assessment/NCLDV_VIRUS_149_hmm", mode: 'symlink'

    input:
    tuple val(id), path(faa)

    output:
    tuple val(id), path("${id}.149_hmm.tblout"), emit:'tblout'

    script:
    """
    hmmsearch --tblout ${id}.149_hmm.tblout --notextw -E 1e-50 --cpu ${task.cpus} ${params.additional_149_NCVOG_hmm} ${faa}
    """
}

/* the ivy's docs said 5 hallmark genes */
/* original code said 8 hallmark genes */
process hmmsearch_with_hallmark_genes {
    publishDir "${params.out}/hallmark/hmm", mode: 'symlink'

    input:
    tuple val(id), path(faa)

    output:
    tuple val(id), path("${id}.hallmark_hmm.tblout"), emit:'tblout'

    script:
    """
    hmmsearch --tblout ${id}.hallmark_hmm.tblout --notextw --cpu ${task.cpus} ${params.hallmark_hmm} ${faa}
    """
}

process detect_hallmark_NCLDV_bin {
    input:
    tuple val(id), path(tblout)

    output:
    tuple val(id), path("${id}.hallmark_presense.tsv"), emit: 'bin'

    script:
    """
    python ${baseDir}/bin/detect_hallmark.py -t ${tblout} -n ${params.hallmark_genes} -s ${params.hallmark_hmm_score_threshold}  -o ${id}.hallmark_presense.tsv
    """
}

process summarize_detected_hallmark_genes {
    publishDir "${params.out}/hallmark", mode: 'symlink'

    input:
    path("table/*")
    
    output:
    path("hallmark_genes.tsv")

    script:
    """
    python ${baseDir}/bin/summarize_table.py -i ./table -p .hallmark_presense.tsv -o hallmark_genes.tsv
    """
}