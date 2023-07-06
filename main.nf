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

/* QC settings */
params.virsorter_groups="dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae"
params.CAT_DB="$baseDir/data/CAT/CAT_prepare_20210107/2021-01-07_CAT_database"
params.CAT_Taxonomy="$baseDir/data/CAT/CAT_prepare_20210107/2021-01-07_taxonomy"
params.additional_NCLDV_149_hmm="$baseDir/data/hmm/NCLDV_149/NCLDV_149.hmm"

params.hallmark_hmm="$baseDir/data/hmm/hallmark/hallmark.hmm"
params.hallmark_genes="DNApolB,MCP_NCLDVs,pATPase_all,Primase_all,RNAP-a_all,RNAP-b_all,TFIIS,VLTF3"
params.hallmark_hmm_score_threshold="150,80,80,80,200,200,100,80"

/* validation settings */
params.hallmark_genes_for_validation="DNApolB,MCP_NCLDVs,pATPase_all,TFIIS,VLTF3"


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

    /*02 detect putative NCLDV bins*/

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

    putative_ncldv_bin_and_prot_ch = classify_NCLDV_bin.out /* filter putative NCLDV bins */
                                                    .filter {
                                                        ncldv: it[2].readLines().last().split('\t').last().toFloat() > params.core_gene_index
                                                    }
                                                    .map{[it[0], it[1]]}.combine(prodigal.out.faa, by: 0)

    classifier_result_list = classify_NCLDV_bin.out.map{it[2]}.toList()

    summarize_NCVOG_results(classifier_result_list)
    /*02 finnished*/

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

    hmmsearch_with_NCLDV_149_hmm(
        putative_ncldv_prot_ch
    )

    assessment_results_ch = putative_ncldv_bin_ch.combine(viralrecall.out.tsv, by: 0)
                                                 .combine(virsorter2.out.tsv, by: 0)
                                                 .combine(CAT.out.txt, by: 0)
                                                 .combine(hmmsearch_with_NCLDV_149_hmm.out.tblout, by: 0)

    summarize_assessment(
        assessment_results_ch
    )

    hmmsearch_with_hallmark_genes(
        putative_ncldv_prot_ch
    )

    detect_hallmark_genes_from_bin(
        hmmsearch_with_hallmark_genes.out.tblout
    )

    hallmark_genes_detection_result_list = detect_hallmark_genes_from_bin.out.table.map{it[1]}.toList()
    
    summarize_detected_hallmark_genes(
        hallmark_genes_detection_result_list
    )
    /*03 finnished*/

    /*04 vaildate NCLDV bins*/
    validate_ncldv_bin_input_ch = summarize_assessment.out.summary.combine(detect_hallmark_genes_from_bin.out.table, by:0)

    validate_ncldv_bin(
        validate_ncldv_bin_input_ch
    )

    validated_ncldv_bin_ch = validate_ncldv_bin.out.table /* filter validated NCLDV bins */
                                                .filter {
                                                    it[1].readLines().last().split('\t').last().toBoolean()
                                                }
                                                .map{[it[0]]}
                                                .combine(putative_ncldv_bin_ch, by: 0)

    validate_ncldv_bin_result_list = validate_ncldv_bin.out.table.map{it[1]}.toList()
    
    summarize_ncldv_bin_validation(
        validate_ncldv_bin_result_list
    )
    /*04 finnished*/

    /*05 remove cellular contig */
    remove_cellular_contig_input_ch = validated_ncldv_bin_ch.combine(summarize_assessment.out.summary, by: 0)

    remove_cellular_contig(
        remove_cellular_contig_input_ch
    )
}

/*
WIP: fastp
WIP: -3 -W 6 -M 30 -q 20 -u 50 -n 0 -p -l 50
*/
process fastp {
    publishDir "${params.out}/binning/quality_check", mode: 'symlink'

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
    publishDir "${params.out}/binning/assembly", mode: 'symlink'

    input:
    tuple val(id), path(forward), path(backward)

    output:
    tuple val(id), path("${id}.megahit.contigs.fa"), emit: 'contig'

    script:
    mem = (task.memory =~ /(\d+).GB/)[0][1] 
    """
    megahit -1 ${forward} -2 ${backward} -m 0.3 --k-min 21 --k-max 141 --k-step 12 -t ${task.cpus} -o megahit --out-prefix ${id}.megahit --min-contig-len 2500
    mv megahit/${id}.megahit.contig.fa .
    """
}

/*
replaced bowtie2 and samtools with coverM, because it's simpler
*/
process coverm {
    publishDir "${params.out}/binning/coverm", mode: 'symlink'
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
    publishDir "${params.out}/binning/metabat2", mode: 'symlink'

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
    publishDir "${params.out}/binning/bins/fasta", mode: 'symlink'

    input:
    path(bin)

    output:
    tuple val(bin_id), path("${bin_id}.fasta"), emit: 'bin'

    script:
    /* define new header */
    run_id=bin.getSimpleName()
    bin_id=bin.getBaseName().replaceAll(/\./,"_")
    """
    seqkit replace -p "k141" -r '${run_id}_k141' ${bin} > ${bin_id}.fasta 
    """
}

/* replace with prodigal-gv? */
process prodigal {
    publishDir "${params.out}/binning/bins/prodigal", mode: 'symlink'

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
    publishDir "${params.out}/detect_NCLDV/tblout", mode: 'symlink'

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
    publishDir "${params.out}/detect_NCLDV/", mode: 'symlink'

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
    publishDir "${params.out}/validate_NCLDV/assessment/viralrecall", mode: 'symlink'

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
    publishDir "${params.out}/validate_NCLDV/assessment/virsorter2", mode: 'symlink'

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
    publishDir "${params.out}/validate_NCLDV/assessment/CAT", mode: 'symlink'

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
process hmmsearch_with_NCLDV_149_hmm {
    publishDir "${params.out}/validate_NCLDV/assessment/NCLDV_149_hmm", mode: 'symlink'

    input:
    tuple val(id), path(faa)

    output:
    tuple val(id), path("${id}.149_hmm.tblout"), emit:'tblout'

    script:
    """
    hmmsearch --tblout ${id}.149_hmm.tblout --notextw -E 1e-50 --cpu ${task.cpus} ${params.additional_NCLDV_149_hmm} ${faa}
    """
}

process summarize_assessment {
    publishDir "${params.out}/validate_NCLDV/assessment/summary", mode: 'symlink'

    input:
    tuple val(id), path(bin), path(viralrecall), path(virsorter2), path(CAT), path(tblout)
    
    output:
    tuple val(id), path("${id}.NCLDV_assessment.tsv"), emit:'summary'

    script:
    """
    python ${baseDir}/bin/summarize_assessment_results.py -f ${bin} -r ${viralrecall} -s ${virsorter2} -c ${CAT} -n ${tblout} -o ${id}.NCLDV_assessment.tsv
    """
}

process hmmsearch_with_hallmark_genes {
    publishDir "${params.out}/validate_NCLDV/hallmark/tblout", mode: 'symlink'

    input:
    tuple val(id), path(faa)

    output:
    tuple val(id), path("${id}.hallmark_hmm.tblout"), emit:'tblout'

    script:
    """
    hmmsearch --tblout ${id}.hallmark_hmm.tblout --notextw --cpu ${task.cpus} ${params.hallmark_hmm} ${faa}
    """
}

process detect_hallmark_genes_from_bin{
    input:
    tuple val(id), path(tblout)

    output:
    tuple val(id), path("${id}.hallmark_presense.tsv"), emit:'table'

    script:
    """
    python ${baseDir}/bin/detect_hallmark.py -t ${tblout} -n ${params.hallmark_genes} -s ${params.hallmark_hmm_score_threshold}  -o ${id}.hallmark_presense.tsv
    """
}

process summarize_detected_hallmark_genes {
    publishDir "${params.out}/validate_NCLDV/hallmark", mode: 'symlink'

    input:
    path("table/*")
    
    output:
    path("hallmark_genes.tsv")

    script:
    """
    python ${baseDir}/bin/summarize_table.py -i ./table -p .hallmark_presense.tsv -o hallmark_genes.tsv
    """
}

process validate_ncldv_bin {
    input:
    tuple val(id), path(assessment_summary), path(hallmark_summary)

    output:
    tuple val(id), path("${id}.NCLDV_validation.tsv"), emit:'table'

    script:
    """
    python ${baseDir}/bin/validate_NCLDV_bin.py -a ${assessment_summary} -m ${hallmark_summary} -c ${params.hallmark_genes_for_validation} -o ${id}.NCLDV_validation.tsv
    """
}

process summarize_ncldv_bin_validation {
    publishDir "${params.out}/validate_NCLDV", mode: 'symlink'

    input:
    path("table/*")
    
    output:
    path("NCLDV_validation.tsv")

    script:
    """
    python ${baseDir}/bin/summarize_table.py -i ./table -p .NCLDV_validation.tsv -o NCLDV_validation.tsv
    """
}

process remove_cellular_contig {
    publishDir "${params.out}/decontamination"

    input:
    tuple val(id), path(bin), path(summary)

    output:
    tuple val(id), path("${id}.decontaminated.fasta")

    script:
    """
    cat ${summary} | awk -F '\t' 'NR>2 && \$6 != 0 {print \$1}' | seqkit grep -f - ${bin} > ${id}.decontaminated.fasta
    """
}