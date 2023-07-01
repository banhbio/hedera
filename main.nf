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

log.info"""
"""


workflow {
    /*01_assembly_and_binning*/
                                
    read_ch = Channel.fromFilePairs("${params.rawread_dir}/*_{1,2}.fastq.gz", flat:true)

    fastp(
        read_ch
    )

    megahit(
        fastp.out.read
    )
    
    forward_reads_list = read_ch.map{it[1]}.toList()
    backward_reads_list = read_ch.map{it[1]}.toList()

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
}

/*
fastp or timommatic?
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
    path("seq1_*.fq.gz")
    path("seq2_*.fq.gz")

    output:
    tuple val(id), path("bam/*"), emit: 'bam'

    script:
    """
    coverm make -r ${contig} -1 seq1_*.fq.gz -2 seq2_*.fq.gz -o ${id} -t ${task.cpus}
    """
}

process metabat2 {
    publishDir "${params.out}/metabat2", mode: 'symlink'

    input:
    tuple val(id), path(contig), path("bam/*")

    output:
    path("${id}.metabat2bin*"), emit: bins
    path("${id}.depth.txt"), emit: depth

    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth ${id}.depth.txt bam/*.bam
    metabat2 -i ${contig} -a ${id}.depth.txt -v -o ${id}.matabat2bin
    """
}
