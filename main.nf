#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.samples_tsv = params.samples_tsv ?: './samples.tsv'
params.db = ""
params.evalue = '0.1'
params.topHits = '5'
params.cpus = '18'
params.publish_dir = './'
params.ref = ''
params.db_name = params.db ? file(params.db).name : ''
params.db_path = params.db ? file(params.db).parent : ''

// Tracy-based consensus assembly process with built-in TSV header check
process tracyAssemble {
    tag "${sampleId}"
    
    conda 'bioconda::tracy'

    publishDir (
        path: "${params.publish_dir}/temp",
        mode: 'copy',
        overwrite: 'true',
    )

    input:
        tuple val(sampleId), path(forward), path(reverse)

    output:
        tuple val(sampleId), path("*.consensus.fasta", optional: true), emit: consensus
        path "*.FAILED", emit: failed

    script:
        """
        if [ -n "${params.ref}" ]; then
            tracy assemble -o ${sampleId} -r ${params.ref} $forward $reverse > "${sampleId}.tracy.log" 2>&1 || true
        else
            tracy assemble -o ${sampleId} $forward $reverse > "${sampleId}.tracy.log" 2>&1 || true
        fi

        # Check for specific error message using case-insensitive regex
        if grep -qiE "at least 2 traces are required" "${sampleId}.tracy.log"; then
            echo "ERROR: Traces don't match for ${sampleId} - Tracy failed" >> "${sampleId}.FAILED"
            exit 0  # Skip this sample gracefully, continue pipeline
        fi

        mv ${sampleId}.cons.fa ${sampleId}.consensus.fasta || exit 1
        """
}

// BLAST process
process blastn {
    tag "${query.baseName}"
    conda 'bioconda::blast'
    
    publishDir "${params.publish_dir}/results", mode: 'copy'

    input:
    tuple val(query_id), path(query)
    path db

    output:
    path ("${query.baseName}_blast_sort.tsv"), emit: blastHits
    path ("${query.baseName}_blast.xml")
    path ("${query.baseName}_blast.html")

    script:
        """
        if [ ! -d "${db}/" ]; then
            exit 1
        fi

        blastn \\
            -query ${query} -db "${db}/${params.db_name}" \\
            -outfmt 11 -out ${query.baseName}_blast.asn \\
            -evalue ${params.evalue} \\
            -num_threads ${params.cpus} || exit 1

        blast_formatter \\
            -archive ${query.baseName}_blast.asn \\
            -outfmt 5 -out ${query.baseName}_blast.xml || exit 1

        blast_formatter \\
            -archive ${query.baseName}_blast.asn \\
            -html -out ${query.baseName}_blast.html || exit 1

        blast_formatter \\
            -archive ${query.baseName}_blast.asn \\
            -outfmt "6 qaccver saccver pident length evalue bitscore mismatch gapopen qstart qend sstart send stitle" -out ${query.baseName}_blast_unsort.tsv || exit 1

        sort -n -r -k 6 ${query.baseName}_blast_unsort.tsv > ${query.baseName}_blast_sort.tsv || exit 1
        """
}

// Top hits process
process blastnTopHits {
    tag "${query.baseName}"
    
    publishDir "${params.publish_dir}/temp", mode: 'copy'

    input:
    path(query)

    output:
    path ("${query.baseName}_topHits.tsv"), emit: topHits

    script:
        """
        if [ ! -s ${query} ]; then
            exit 1
        fi

        head -n ${params.topHits} ${query} > ${query.baseName}_topHits.tsv || exit 1
        """
}

// Concatenate hits process
process blastHits_concat {

    publishDir "${params.publish_dir}/results", mode: 'copy'

    input:
    file combinedFiles

    output:
    path "concatenatedHits_withHeaders.tsv"

    script:

        """
        if [ ${combinedFiles.size()} -eq 0 ]; then
            exit 1
        fi

        cat ${combinedFiles} >> concatenatedHits.tsv || exit 1
        awk 'BEGIN{print "qaccver\\tsaccver\\tpident\\tlength\\tevalue\\tbitscore\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tstitle"}1' concatenatedHits.tsv > concatenatedHits_withHeaders.tsv || exit 1
        """
}

// Workflow
workflow {
    log.info """

    s a n g e r F L O W
    A Sanger sequencing-based bioinformatics pipeline for pest and pathogen diagnosis

    Asad Prodhan
    =======================================

    """
    .stripIndent()

    // Read samples.tsv: sampleId, read1 (forward SCF), read2 (reverse SCF)
    // First, check the header manually before parsing
    def header = new File("${params.samples_tsv}").readLines().getAt(0)
    if (header != "sampleId\tread1\tread2") {
        throw new IllegalArgumentException("ERROR: samples.tsv header must be exactly: sampleId\\tread1\\tread2\nCurrent header: '$header'")
    }

    def samplesChannel = channel.fromPath(params.samples_tsv)
        | splitCsv( header: true, sep: '\t' )
        | map { row -> tuple( row.sampleId, file(row.read1), file(row.read2) ) }

    // Run Tracy assembly
    tracyAssemble( samplesChannel )

    // Run BLAST pipeline (requires --db parameter) only on successful assemblies
    blastn( tracyAssemble.out.consensus, params.db_path )
    blastnTopHits(blastn.out.blastHits)
    blastHits_concat(blastnTopHits.out.topHits.collect())

    // Log failed samples
    tracyAssemble.failed.view { id -> log.error "FAILED: $id" }
}
