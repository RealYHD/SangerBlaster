#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.samples_tsv = params.samples_tsv ?: './samples.tsv'
params.evalue = '0.1'
params.topHits = '5'
params.cpus = '18'
params.publish_dir = './'
params.ref = ''

// Tracy-based consensus assembly process with built-in TSV header check
process tracyAssemble {
    tag "${sampleId}"
    
    conda 'bioconda::tracy'
    // Tracy failing on one sample shouldn't take down the whole run: keep the
    // pipeline alive and let the optional outputs decide what flows downstream.
    errorStrategy 'ignore'

    publishDir (
        path: "${params.publish_dir}/temp",
        mode: 'copy',
        overwrite: 'true',
    )

    input:
        tuple val(sampleId), path(forward), path(reverse)

    output:
        tuple val(sampleId), path("*.consensus.fasta"), optional: true, emit: consensus
        path "*.FAILED", optional: true, emit: failed

    script:
        """
        # Run Tracy; capture its exit code without letting -e abort the script (|| keeps us alive)
        if [ -n "${params.ref}" ]; then
            tracy assemble -o ${sampleId} -r ${params.ref} $forward $reverse > "${sampleId}.tracy.log" 2>&1 || tracy_exit=\$?
        else
            tracy assemble -o ${sampleId} $forward $reverse > "${sampleId}.tracy.log" 2>&1 || tracy_exit=\$?
        fi

        # Move the consensus into place (|| true tolerates a stale/same-inode target)
        if [ -f ${sampleId}.cons.fa ]; then
            mv ${sampleId}.cons.fa ${sampleId}.consensus.fasta || true
        fi

        # Exit with Tracy's real exit code so Nextflow knows this sample failed
        # (errorStrategy 'ignore' still lets the rest of the run continue)
        if [ -n "\${tracy_exit:-}" ]; then
            exit \${tracy_exit}
        fi
        """
}

// BLAST process
// Missing --db index is a config error checked upfront in the workflow; a per-sample
// blastn failure is caught by errorStrategy 'ignore' so the sample is skipped downstream.
process blastn {
    tag "${query.baseName}"
    conda 'bioconda::blast'
    errorStrategy 'ignore'
    
    publishDir "${params.publish_dir}/results", mode: 'copy'

    input:
    tuple val(query_id), path(query)

    output:
    path ("${query.baseName}_blast_sort.tsv"), optional: true, emit: blastHits
    path ("${query.baseName}_blast.xml"), optional: true
    path ("${query.baseName}_blast.html"), optional: true

    script:
        """
        blastn -query ${query} -db "${params.db}" -outfmt 11 -out ${query.baseName}_blast.asn -evalue ${params.evalue} -num_threads ${params.cpus} || true

        blast_formatter -archive ${query.baseName}_blast.asn -outfmt 5 -out ${query.baseName}_blast.xml || true

        blast_formatter -archive ${query.baseName}_blast.asn -html -out ${query.baseName}_blast.html || true

        blast_formatter -archive ${query.baseName}_blast.asn -outfmt "6 qaccver saccver pident length evalue bitscore mismatch gapopen qstart qend sstart send stitle" -out ${query.baseName}_blast_unsort.tsv || true

        # Sort by bitscore
        sort -n -r -k 6 ${query.baseName}_blast_unsort.tsv > ${query.baseName}_blast_sort.tsv || true
        """
}

// Top hits process
// An empty query is a per-sample condition (empty consensus → no hits); caught by
// errorStrategy 'ignore' so downstream simply skips that sample.
process blastnTopHits {
    tag "${query.baseName}"
    errorStrategy 'ignore'
    
    publishDir "${params.publish_dir}/temp", mode: 'copy'


    input:
    path(query)

    output:
    path ("${query.baseName}_topHits.tsv"), emit: topHits

    script:
        """
        head -n ${params.topHits} ${query} > ${query.baseName}_topHits.tsv || true
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
    // Validate required parameters upfront (fail fast)
    if (!params.db) {
        throw new IllegalArgumentException("ERROR: --db parameter is required for BLAST workflow.\nUsage: nextflow run main.nf --db /path/to/database/base\n")
    }

    // Validate the BLAST database index exists (a missing index is a config error, not
    // a per-sample condition, so it stops the run before any process starts)
    def dbIndex = ["${params.db}.ndb", "${params.db}.nhd"]
    if (!dbIndex.any { new File(it).exists() }) {
        throw new IllegalArgumentException("ERROR: BLAST database not found: ${params.db}.ndb or ${params.db}.nhd\n")
    }

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
    blastn( tracyAssemble.out.consensus )
    blastnTopHits(blastn.out.blastHits)
    blastHits_concat(blastnTopHits.out.topHits.collect())
}
