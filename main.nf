#!/usr/bin/env nextflow

/*
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl=2

version = '0.0.5'

log.info """
hedera: This is a pipleline to assemble and binning NCLDV MAGs from metagenomic data.
version:${version}
"""

if (params.help) {
    log.info """
    /* help messeage */
    """
    exit 0
}

/* utils */
def createTable(String items1, String items2, String column1Name, String column2Name) {
    // Split the input strings into lists
    def itemList1 = items1.split(",")
    def itemList2 = items2.split(",")

    // Check if the sizes of the lists match
    if(itemList1.size() != itemList2.size()) {
        throw new IllegalArgumentException("The size of both item lists must be the same.")
    }

    // Find the longest length of the strings in the lists
    def maxColumn1Length = itemList1.collect { it.length() }.max().compareTo(column1Name.length()) > 0 ? itemList1.collect { it.length() }.max() : column1Name.length()
    def maxColumn2Length = itemList2.collect { it.length() }.max().compareTo(column2Name.length()) > 0 ? itemList2.collect { it.length() }.max() : column2Name.length()

    // Make the table header
    String table = "| ${column1Name.padRight(maxColumn1Length)} | ${column2Name.padRight(maxColumn2Length)} |\n" 
    table += "| ${"-".multiply(maxColumn1Length)} | ${"-".multiply(maxColumn2Length)} |\n"

    // Create each row
    for(int i = 0; i < itemList1.size(); i++) {
        table += "| ${itemList1[i].padRight(maxColumn1Length)} | ${itemList2[i].padRight(maxColumn2Length)} |\n"
    }

    return table
}

NCVOG_table = createTable(params.conserved_20_NCVOGs, params.conserved_20_NCVOG_weights, "NCVOG", "Weight")
hallmark_gene_table = createTable(params.hallmark_genes, params.hallmark_hmm_score_threshold, "Hallmark gene", "Score threshold")

/* Print Infomation */
log.info"""
General parameters in this run:
"""
if(!params.from_bins){
    if(params.after_qc){
        log.info"""
            --after_qc : true
                Input reads are treated as quality-checked (e.g., with fastp, trimmomatic, etc.).
        """
    }

    if(params.from_contig){
        log.info"""
            --from_contig : true
                The assembly step will skipped. Please ensure --input_contigs is not left empty.
        """
        if(params.input_contigs.isEmpty()){
                log.info""
                exit 1, "Oops! --input_contigs is empty."
        }
    }

    log.info"""
        Input reads  : ${params.input_reads}
    """

    if(!params.after_qc){
        log.info"""
            Quality check settings:
                fastp parameter : ${params.fastp_parameter} 
        """
    }

    if(params.from_contig){
        log.info"""
            Input contigs: ${params.input_contigs}
        """
    }else{
        log.info"""
            Assembly settings:
                Assembly name             : ${params.assembly_name}
                Memory per machine' total : ${params.assembly_memory_per_machine}
                Kmer parameter            : ${params.assembly_kmer}
                Minimun contig length     : ${params.assembly_min_contig_len}
        """
    }
    log.info"""
        Binnig settings:
            Minimum contig length resucued from binning leftovers: ${params.leftover_length}
    """
}else{
    log.info"""
        --from_bins : true
            Start the pipleline from bins. Please ensure --input_bins and --input_depth are not left empty.
        
        Input bins  : ${params.input_bins}
        Input depth : ${params.input_depth}
   """
}
    

log.info"""
    Detect NCLDV bin settings:
        Core gene (20 NCVOG) hmms : ${params.conserved_20_NCVOG_hmm}
        Core gene index           : ${params.core_gene_index}
        Weights of NCVOGs         :

${NCVOG_table.split('\n').collect { '       ' + it }.join('\n')}
    
    Assessment NCLDV bin settings:
        Virsoter groups             : ${params.virsorter_groups}
        CAT DB                      : ${params.CAT_DB}
        CAT Taxonomy                : ${params.CAT_Taxonomy}
        Additional 149 hmms         : ${params.additional_NCLDV_149_hmm} 
        Additional 149 hmms e-value : ${params.additional_NCLDV_149_hmm_evalue} 
        NCLDV hallmark genes        :

${hallmark_gene_table.split('\n').collect { '       ' + it }.join('\n')}

    Validate NCLDV bin settings:
        NCLDV hallmark genes : ${params.hallmark_genes_for_validation}

    Delineage NCLDV bin settings:
        NCLDV single copy genes : ${params.hallmark_scgs}
"""



workflow {
    /*01 quality control of reads*/

    if(!params.from_bins) {
        read_ch = Channel.fromFilePairs(params.input_reads, flat:true)
        if(!params.after_qc){
            fastp(
                read_ch
            )
            qc_read_ch = fastp.out.read
        }else{
            qc_read_ch = read_ch
        }
    /*01 finnish*/

    /*02 assembly*/
        forward_reads_list = qc_read_ch.map{it[1]}.toList()
        backward_reads_list = qc_read_ch.map{it[2]}.toList()

        if(!params.from_contig){
            megahit(
                forward_reads_list,
                backward_reads_list
            )

            contig_ch = megahit.out.contig
        }else{
            contig_ch = Channel.fromPath(params.input_contigs)
        }
    /*02 finnish*/

    /*03 binnning*/
        coverm(
            contig_ch,
            forward_reads_list,
            backward_reads_list
        )

        summarize_depth(
            coverm.out.bam
        )

        depth_ch = summarize_depth.out.depth

        metabat2(
            contig_ch,
            depth_ch
        )

        bin_ch = metabat2.out.bins.flatten().mix(metabat2.out.single_contig.flatten())
    }else{
        bin_ch = Channel.fromPath(params.input_bins)
        depth_ch = Channel.fromPath(params.input_depth)
    }

    rename_bin(
        bin_ch
    )
    /*03 finnish*/

    /*04 filter putative NCLDV bins*/

    /* collect all bins */
    prodigal(
        rename_bin.out.bin
    )

    hmmsearch_with_NCVOGs(
        prodigal.out.faa
    )

    classifier_input_ch = rename_bin.out.bin.combine(hmmsearch_with_NCVOGs.out.tblout, by: 0)
    
    classify_NCLDV_bin(
        classifier_input_ch
    )

    putative_ncldv_bin_and_prot_ch = classify_NCLDV_bin.out /* filter putative NCLDV bins */
                                                    .filter {
                                                        ncldv: it[2].readLines().last().split('\t').last().toFloat() > params.core_gene_index
                                                    }
                                                    .map{[it[0], it[1]]}.combine(prodigal.out.faa, by: 0)
    /*04 finnish*/

    /*05 validate putative NCLDV bins*/
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

    detect_hallmark_genes_from_bin_ch = putative_ncldv_bin_ch.combine(hmmsearch_with_hallmark_genes.out.tblout, by: 0)

    detect_hallmark_genes_from_bin(
        detect_hallmark_genes_from_bin_ch
    )

    validate_ncldv_bin_input_ch = summarize_assessment.out.summary.combine(detect_hallmark_genes_from_bin.out.table, by: 0)

    validate_ncldv_bin(
        validate_ncldv_bin_input_ch
    )

    validated_ncldv_bin_ch = validate_ncldv_bin.out.table /* filter validated NCLDV bins */
                                                .filter {
                                                    validated: it[1].readLines().last().split('\t').last().toBoolean()
                                                }
                                                .combine(putative_ncldv_bin_ch, by: 0)
                                                .map{[it[0], it[2]]}
    /*05 finnish*/

    /*06 decontaminate cellular contig */
    remove_cellular_contig_input_ch = validated_ncldv_bin_ch.combine(summarize_assessment.out.summary, by: 0)

    remove_cellular_contig(
        remove_cellular_contig_input_ch
    )
    /*06 finnish*/

    /*07 delineage NCLDV bin*/

    delineage_candidate_input_ch = remove_cellular_contig.out.bin
                                                        .combine(detect_hallmark_genes_from_bin.out.table_per_contig, by: 0)
                                                        .combine(depth_ch)

    find_delineage_candidate(
        delineage_candidate_input_ch,
    )

    count_tetramer(
        remove_cellular_contig.out.bin
    )

    find_delineage_candidate.out.candidate
                            .combine(count_tetramer.out.tetramer, by:0)
                            .branch{
                                    candidate: it[3].readLines().last().split('\t').last().toBoolean()
                                    clean: true
                                }
                            .set {delineage}

    delineage_clean_bin_ch = delineage.clean.map{[it[0], it[1], it[2], it[4]]}
    
    /*
        Delineage process requires 4 files: 1) bin fasta.
                                           2) depth of contigs in bin (depth.txt).
                                           3) tetra nucleotide frequency.
                                           3) NCLDV score of contigs in bin (summary.tsv).
                                           4) Taxonomic annotation of contigs in bin (CAT.txt).
    */
    delineage_input_ch = delineage.candidate.map{[it[0], it[1], it[2], it[4]]}
                                        .combine(detect_hallmark_genes_from_bin.out.table_per_contig, by: 0)
                                        .combine(summarize_assessment.out.summary, by: 0)

    delineage_bin(
        delineage_input_ch
    )

    postdelineage_input_ch = delineage_bin.out.bin.flatMap{
                                                    item ->
                                                        /* check if the delineaged output is a list */
                                                        if (item[0] instanceof List) {
                                                            item[0].collect {bin -> [bin, item.tail()].flatten()}
                                                        }else{
                                                            // If the bin is single, return the item
                                                            [item]
                                                        }   
                                                    }

    postdelineage(
        postdelineage_input_ch
    )

    /*07 finnish */

    /*08 seconde decontamination*/
    after_delineage_bin_ch = delineage_clean_bin_ch.combine(summarize_assessment.out.summary, by:0)
                                                   .mix(postdelineage.out.bin)

    second_decontamination(
        after_delineage_bin_ch
    )

    /*08 finish*/

    /*09 summary table */
    table_ch = classify_NCLDV_bin.out.map{it[2]}.toList().map {["04_filter_NCLDV", "20NCVOG_weight.tsv" , it]}
                .mix(detect_hallmark_genes_from_bin.out.table.map{it[1]}.toList().map {["05_validate_NCLDV/hallmark", "hallmark_gene.tsv", it]})
                .mix(validate_ncldv_bin.out.table.map{it[1]}.toList().map {["05_validate_NCLDV", "NCLDV_validation.tsv", it]})
                .mix(find_delineage_candidate.out.candidate.map{it[3]}.toList().map {["07_delineage_NCLDV/", "delineage_candidate.tsv", it]})
                .mix(delineage_bin.out.table.toList().map {["07_delineage_NCLDV/delineage","delineage_summary.tsv", it]})
    
    summarize_table(
        table_ch
    )
    /*09 finnish*/
}

process fastp {
    publishDir "${params.out}/01_qc_reads/fastq", pattern: '*.fastq.gz', mode: 'symlink'
    publishDir "${params.out}/01_qc_reads/html", pattern: '*.html', mode: 'symlink'
    publishDir "${params.out}/01_qc_reads/json", pattern: '*.json', mode: 'symlink'

    input:
    tuple val(id), path("seq1.fastq.gz"), path("seq2.fastq.gz")

    output:
    tuple val("${id}"), path("${id}_qcd_1.fastq.gz"), path("${id}_qcd_2.fastq.gz"), emit: 'read'
    path("${id}_fastp_report.html")
    path("${id}_fastp_report.json")

    script:
    """
    fastp -i seq1.fastq.gz -I seq2.fastq.gz -o ${id}_qcd_1.fastq.gz -O ${id}_qcd_2.fastq.gz -w ${task.cpus} -h ${id}_fastp_report.html -j ${id}_fastp_report.json ${params.fastp_parameter}
    """
}

/*
rename contig header by run id
*/
process megahit {
    publishDir "${params.out}/02_assembly", mode: 'symlink'

    input:
    path("seq*_1.fastq.gz")    
    path("seq*_2.fastq.gz")    

    output:
    path("${id}.megahit.rename.contigs.fa"), emit: 'contig'

    script:
    id = params.assembly_name
    mem = (task.memory =~ /(\d+).GB/)[0][1] 
    """
    forward=\$(ls seq*_1.fastq.gz | tr '\n' ',' | sed 's/,\$//')
    backward=\$(ls seq*_2.fastq.gz | tr '\n' ',' | sed 's/,\$//')
    megahit -1 \$forward -2 \$backward -m ${params.assembly_memory_per_machine}  -t ${task.cpus} -o megahit --out-prefix ${id}.megahit ${params.assembly_kmer} --min-contig-len ${params.assembly_min_contig_len}
    seqkit replace -p '^' -r '${id}_' megahit/${id}.megahit.contigs.fa > ${id}.megahit.rename.contigs.fa
    """
}

process coverm {
    publishDir "${params.out}/03_binning/coverm", mode: 'symlink'

    input:
    path(contig)
    path("seq1/*")
    path("seq2/*")

    output:
    path("bam/*"), emit: 'bam'

    script:
    """
    coverm make -r ${contig} -1 seq1/* -2 seq2/* -o bam -t ${task.cpus}
    """
}

process summarize_depth {
    publishDir "${params.out}/03_binning/depth", mode: 'symlink'

    input:
    path("bam/*")

    output:
    path("depth.txt"), emit: 'depth'

    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth depth.txt bam/*.bam
    """
}

process metabat2 {
    publishDir "${params.out}/03_binning/", mode: 'symlink'

    input:
    path(contig)
    path(depth)

    output:
    path("metabat2bin/*"), emit: 'bins'
    path("metabat2sc/*"), emit: 'single_contig', optional: true

    script:
    """
    metabat2 -i ${contig} -a ${depth} -t ${task.cpus} -v -o metabat2bin/metabat2bin
    cat metabat2bin/* | seqkit seq -ni | seqkit grep -v -f - ${contig} | seqkit replace -p "\s.+" | seqkit seq -m ${params.leftover_length} > metabat2sc.fasta
    mkdir metabat2sc
    [ -s "metabat2sc.fasta" ] || exit 0
    seqkit split -s 1 -O tmp metabat2sc.fasta
    for old in tmp/*; do base=\$(basename \$old); suffix=\${base#*.}; number_fasta=\$(echo \$suffix | cut -d'_' -f 2 | sed 's/^0*//'); mv "\$old" metabat2sc/metabat2sc."\${number_fasta}"; done
    rm -rf tmp
    """
}


/* This is to avoid viralrecall error*/
process rename_bin {
    publishDir "${params.out}/03_binning/bin/fasta", mode: 'symlink'

    input:
    path(bin)

    output:
    tuple val(id), path("${id}.fasta"), emit: 'bin'

    script:
    id = bin.getBaseName().replace(".", "_")
    """
    cp ${bin} ${id}.fasta || exit 0
    """
}

process prodigal {
    publishDir "${params.out}/03_binning/bin/prodigal", mode: 'symlink'

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
    publishDir "${params.out}/04_filter_NCLDV/tblout", mode: 'symlink'

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
    python ${baseDir}/bin/classify_bin_20NCVOG.py -b ${id} -f ${bin} -t ${tblout} -n ${params.conserved_20_NCVOGs} -w ${params.conserved_20_NCVOG_weights} -o ${id}.20NCVOG_weight.tsv 
    """
}

process viralrecall {
    publishDir "${params.out}/05_validate_NCLDV/assessment/viralrecall", mode: 'symlink'

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
    publishDir "${params.out}/05_validate_NCLDV/assessment/virsorter2", mode: 'symlink'

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
    publishDir "${params.out}/05_validate_NCLDV/assessment/CAT", mode: 'symlink'

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

process hmmsearch_with_NCLDV_149_hmm {
    publishDir "${params.out}/05_validate_NCLDV/assessment/NCLDV_149_hmm", mode: 'symlink'

    input:
    tuple val(id), path(faa)

    output:
    tuple val(id), path("${id}.149_hmm.tblout"), emit:'tblout'

    script:
    """
    hmmsearch --tblout ${id}.149_hmm.tblout --notextw -E ${params.additional_NCLDV_149_hmm_evalue} --cpu ${task.cpus} ${params.additional_NCLDV_149_hmm} ${faa}
    """
}

process summarize_assessment {
    publishDir "${params.out}/05_validate_NCLDV/summary", mode: 'symlink'

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
    publishDir "${params.out}/05_validate_NCLDV/hallmark/tblout", mode: 'symlink'

    input:
    tuple val(id), path(faa)

    output:
    tuple val(id), path("${id}.hallmark_hmm.tblout"), emit:'tblout'

    script:
    """
    hmmsearch --tblout ${id}.hallmark_hmm.tblout --notextw --cpu ${task.cpus} ${params.hallmark_hmm} ${faa}
    """
}

process detect_hallmark_genes_from_bin {
    input:
    tuple val(id), path(bin), path(tblout)

    output:
    tuple val(id), path("${id}.hallmark_total.tsv"), emit:'table'
    tuple val(id), path("${id}.hallmark_per_contigs.tsv"), emit:'table_per_contig'

    script:
    """
    python ${baseDir}/bin/detect_hallmark.py -b ${id} -f ${bin} -t ${tblout} -n ${params.hallmark_genes} -s ${params.hallmark_hmm_score_threshold}  -o ${id}.hallmark_total.tsv -O ${id}.hallmark_per_contigs.tsv
    """
}

process validate_ncldv_bin {
    input:
    tuple val(id), path(assessment_summary), path(hallmark_summary)

    output:
    tuple val(id), path("${id}.NCLDV_validation.tsv"), emit:'table'

    script:
    """
    python ${baseDir}/bin/validate_NCLDV_bin.py -b ${id} -a ${assessment_summary} -m ${hallmark_summary} -c ${params.hallmark_genes_for_validation} -o ${id}.NCLDV_validation.tsv
    """
}

process remove_cellular_contig {
    publishDir "${params.out}/06_decontaminate_NCLDV/bin", mode: 'symlink'

    input:
    tuple val(id), path(bin), path(summary)

    output:
    tuple val(id), path("${id}.decontaminated.fasta"), emit:'bin'

    script:
    """
    cat ${summary} | csvtk filter -t -f "NCLDV_score>0" | sed '1d' | cut -f1 | seqkit grep -f - ${bin} > ${id}.decontaminated.fasta
    """
}

process find_delineage_candidate {
    input:
    tuple val(id), path(bin), path(hallmark_summary), path(depth)

    output:
    tuple val(id), path(bin), path("${id}.depth.txt"), path("${id}.delineage_candidate.tsv"), emit:'candidate'

    script:
    """
    cat ${bin} | seqkit seq -ni | csvtk grep -t -f1 -P - ${depth} > ${id}.depth.txt
    python ${baseDir}/bin/find_delineage_candidate.py -b ${id} -m ${hallmark_summary} -s ${params.hallmark_scgs} -d ${id}.depth.txt -o ${id}.delineage_candidate.tsv
    """
}

process count_tetramer {
    input:
    tuple val(id), path(bin)

    output:
    tuple val(id), path("${id}.tetramer.tsv"), emit:'tetramer'

    script:
    """
    cat ${bin} | seqkit replace -p "\s.+" | cgat fasta2kmercontent -k 4 -p | grep -v "^#" > ${id}.tetramer.tsv
    """
}

process delineage_bin {
    publishDir "${params.out}/07_delineage_NCLDV/delineage", pattern: "bin/*", mode: 'symlink'
    publishDir "${params.out}/07_delineage_NCLDV/delineage/tree", pattern: "*.tree", mode: 'symlink'

    input:
    tuple val(id), path(bin), path(depth), path(tetramer), path(hallmark_summary), path(assessment)

    output:
    tuple path("bin/*"), path(depth), path(assessment),  emit:'bin'
    path("${id}.tree"), emit:'tree'
    path("${id}.delineage_summary.tsv"), emit: 'table'

    script:
    """
    mkdir -p bin
    python ${baseDir}/bin/delineage.py -b ${id} \
                                       -f ${bin} \
                                       -t ${tetramer} \
                                       -d ${depth} \
                                       -m ${hallmark_summary} \
                                       -n ${assessment} \
                                       -s ${params.hallmark_scgs} \
                                       -o bin \
                                       -O ${id}.tree \
                                       -S ${id}.delineage_summary.tsv
    """
}

process postdelineage {
    input:
    tuple path(bin), path(depth), path(assessment)

    output:
    tuple val(id), path(bin), path("${id}.depth.txt"), path("${id}.tetramer.tsv"), path("${id}.NCLDV_assessment.txt"), emit:'bin'

    script:
    id = bin.getSimpleName()
    """
    cat ${bin} | seqkit replace -p "\s.+" | cgat fasta2kmercontent -k 4 -p | grep -v "^#" > ${id}.tetramer.tsv
    cat ${bin} | seqkit seq -ni | csvtk grep -t -f1 -P - ${depth} > ${id}.depth.txt
    cat ${bin} | seqkit seq -ni | csvtk grep -t -f1 -P - ${assessment} > ${id}.NCLDV_assessment.txt
    """
}

process second_decontamination {
    publishDir "${params.out}/09_final_NCLDV_MAG", pattern: "*.fasta", mode: 'symlink'
    publishDir "${params.out}/08_decontaminate_NCLDV_2nd/summary", pattern: "*.tsv",  mode: 'symlink'

    input:
    tuple val(id), path(bin), path(depth), path(tetramer), path(assessment_summary)

    output:
    tuple val(id), path(bin), emit:'MAG'
    path("${id}.second_decontamination.tsv"), emit:'table'

    script:
    """
    python ${baseDir}/bin/second_decontamination.py -f ${bin} -c ${depth} -t ${tetramer} -a ${assessment_summary} -d ${id}.second_decontamination.tsv -o ${id}.mag.fasta
    """
}

process summarize_table {
    publishDir "${params.out}/${path}/", mode: 'symlink'

    input:
    tuple val(path), val(file), path("table/*")
    
    output:
    path("${file}")

    script:
    """
    touch ${file}
    csvtk concat -t table/* > ${file} || exit 0
    """
}
