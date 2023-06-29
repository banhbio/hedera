#https://github.com/OpenGene/fastp
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz

#nextflow script by ban-kun
process fastp {
    publishDir "${params.out}/qc", mode: 'symlink'

    input:
    tuple val(id), path("seq?.fq.gz")

    output:
    tuple val("$id"), path("${id}_forward_qcd.fq.gz"), path("${id}_backward_qcd.fq.gz"), emit: 'read'
    path("${id}_fastp_report.html")
    path("${id}_fastp_report.json")

    script:
    """
    fastp -i seq1.fq.gz -I seq2.fq.gz -o ${id}_forward_qcd.fq.gz -O ${id}_backward_qcd.fq.gz -3 -W 6 -M 30 -q 20 -u 50 -n 0 -p -l 50 -w ${
task.cpus} -h ${id}_fastp_report.html -j ${id}_fastp_report.json
    """
}