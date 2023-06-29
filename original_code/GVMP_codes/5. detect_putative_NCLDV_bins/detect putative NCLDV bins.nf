#workflow: 4 steps
#input: *.fa (bin fast files)
#output: a list of putative NCLDV bins in text

#!/usr/bin/env nextflow
import nextflow.splitter.CsvSplitter
/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl=2

version = '0.1'

params.cpu = 8

params.help = false
params.resume = false

params.inputdir = "$baseDir/test"
params.out="$baseDir/out"

workflow {

    #ids = fetchRunAccessionsAndIDs(params.filereport)
    read_ch = Channel.fromSRA(ids, apiKey: params.api, protocol:'https')

    change_contig_id(
        read_ch
    )

    prodigal(
        fastp.out.read,
        params.phylofrashDB
    )

    hmmsearch(
        fastp.out.read,
        params.kaijuDB
    )

    input_ch = Channel.fromPath("${params.inputdir}/*.fasta")

    NCLDVbinfinder(input_ch)
}

process change_contig_id {
    publishDir "${params.out}/qc", mode: 'symlink'

    input:
    tuple val(id), path("seq?.fq.gz")

    output:
    tuple val("$id"), path("${id}_forward_qcd.fq.gz"), path("${id}_backward_qcd.fq.gz"), emit: 'read'
    path("${id}_fastp_report.html")

    script:
    """
    #contig id 前加上bin id
for f in ./*.fa;
    do
    id=`basename ${f%.*}`
    sed -i "s/k141/$id\_k141/" $f
    done
    #rename metabat "$sid"_metabat $dir*.fa
    """
}