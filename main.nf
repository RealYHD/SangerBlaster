#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.samples_tsv = './samples.tsv'
params.db = ""
params.evalue = '0.1'
params.topHits = '5'
params.cpus = '18'
params.publish_dir = './'
params.ref = ''
params.db_name = params.db ? file(params.db).name : ''
params.db_path = params.db ? file(params.db).parent : ''

// =============================================================================
// Tracy-based consensus: takes forward/reverse SCF chromatograms + reference,
// produces a consensus FASTA. Replaces the old bioawk/Clustal/EMBOSS pipeline.
// =============================================================================
process tracyAssemble {
    tag "${sampleId}"

    publishDir (
        path: "${params.publish_dir}/temp",
        mode: 'copy',
        overwrite: 'true',
    )

    input:
        tuple val(sampleId), path(forward), path(reverse)

    output:
        tuple val(sampleId), path("*.fasta"), emit: consensus

    script:
        """
        # Ensure Tracy binary is available (static binary, downloaded once to workdir)
        if [ ! -f tracy ]; then
            apt-get update -qq && apt-get install -y -qq wget && \
            wget -q -O tracy https://github.com/gear-genomics/tracy/releases/download/v0.9.1/tracy-v0.9.1-linux-amd64 && \
            chmod +x tracy
        fi

        # Reference-guided or de-novo assembly of forward + reverse chromatograms
        # Tracy writes consensus to {outprefix}.cons.fa, not stdout
        if [ -n "${params.ref}" ]; then
            ./tracy assemble -o "${sampleId}" -r ${params.ref} $forward $reverse 2>/dev/null
        else
            ./tracy assemble -o "${sampleId}" $forward $reverse 2>/dev/null
        fi

        # Rename Tracy's output to match Nextflow's *.fasta output pattern
        mv "${sampleId}.cons.fa" "${sampleId}.consensus.fasta"
        """
}

// =============================================================================
// BLAST processes (unchanged from original sangerFlow)
// =============================================================================
process blastn {
    tag "${query.baseName}"
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
        blastn \\
            -query ${query} -db "${db}/${params.db_name}" \\
            -outfmt 11 -out ${query.baseName}_blast.asn \\
            -evalue ${params.evalue} \\
            -num_threads ${params.cpus}

        blast_formatter \\
            -archive ${query.baseName}_blast.asn \\
            -outfmt 5 -out ${query.baseName}_blast.xml

        blast_formatter \\
            -archive ${query.baseName}_blast.asn \\
            -html -out ${query.baseName}_blast.html

        blast_formatter \\
            -archive ${query.baseName}_blast.asn \\
            -outfmt "6 qaccver saccver pident length evalue bitscore mismatch gapopen qstart qend sstart send stitle" -out ${query.baseName}_blast_unsort.tsv

        sort -n -r -k 6 ${query.baseName}_blast_unsort.tsv > ${query.baseName}_blast_sort.tsv
        """
}

process blastnTopHits {
    tag "${query.baseName}"
    publishDir "${params.publish_dir}/temp", mode: 'copy'

    input:
    path(query)

    output:
    path ("${query.baseName}_topHits.tsv"), emit: topHits

    script:
        """
        head -n ${params.topHits} ${query} > ${query.baseName}_topHits.tsv
        """
}

process blastHits_concat {

    publishDir "${params.publish_dir}/results", mode: 'copy'

    input:
    file combinedFiles

    output:
    path "concatenatedHits_withHeaders.tsv"

    script:

        """
        cat ${combinedFiles} >> concatenatedHits.tsv
        awk 'BEGIN{print "qaccver\\tsaccver\\tpident\\tlength\\tevalue\\tbitscore\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tstitle"}1' concatenatedHits.tsv > concatenatedHits_withHeaders.tsv
        """
}

// =============================================================================
// Workflow
// =============================================================================
workflow {
    log.info """\\

    s a n g e r F L O W
    A Sanger sequencing-based bioinformatics pipeline for pest and pathogen diagnosis 

    Asad Prodhan
    =======================================

    """
    .stripIndent()

    // Read samples.tsv: sampleId, read1 (forward SCF), read2 (reverse SCF)
    // Then use Tracy to assemble a consensus FASTA directly from chromatograms
    def samplesChannel = channel.fromPath( params.samples_tsv )
        | splitCsv( header: true, sep: '\t' )
        | map { row -> tuple( row.sampleId, file(row.read1), file(row.read2) ) }

    tracyAssemble( samplesChannel )

    blastn( tracyAssemble.out.consensus, params.db_path )
    blastnTopHits(blastn.out.blastHits)

    blastHits_concat(blastnTopHits.out.topHits.collect())
}
