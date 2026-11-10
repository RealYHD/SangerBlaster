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

// Tracy-based consensus: takes forward/reverse SCF chromatograms + reference,
// produces a consensus FASTA. Replaces the old bioawk/Clustal/EMBOSS pipeline.
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
        # Install wget if not available with timeout
        if ! command -v wget &> /dev/null; then
            timeout 30 apt-get update -qq && timeout 60 apt-get install -y -qq wget || exit 1
        fi

        # Ensure Tracy binary is available (static binary, downloaded once to workdir)
        if [ ! -f tracy ]; then
            # Use timeout with wget or curl, fail after 30 seconds total
            if command -v wget &> /dev/null; then
                timeout 30 wget -q -O tracy https://github.com/gear-genomics/tracy/releases/download/v0.9.1/tracy-v0.9.1-linux-amd64 || exit 1
                chmod +x tracy
            elif command -v curl &> /dev/null; then
                timeout 30 curl -L -o tracy https://github.com/gear-genomics/tracy/releases/download/v0.9.1/tracy-v0.9.1-linux-amd64 || exit 1
                chmod +x tracy
            else
                echo "ERROR: Neither wget nor curl is available to download Tracy!"
                exit 1
            fi
            
            # Verify download succeeded and file is executable
            if [ ! -s tracy ] || [ ! -x tracy ]; then
                echo "ERROR: Failed to download or make Tracy executable!"
                exit 1
            fi
        fi

        # Reference-guided or de-novo assembly of forward + reverse chromatograms
        # Tracy writes consensus to {outprefix}.cons.fa, not stdout
        if [ -n "${params.ref}" ]; then
            timeout 3600 ./tracy assemble -o "${sampleId}" -r ${params.ref} $forward $reverse 2>/dev/null || exit 1
        else
            timeout 3600 ./tracy assemble -o "${sampleId}" $forward $reverse 2>/dev/null || exit 1
        fi

        # Rename Tracy's output to match Nextflow's *.fasta output pattern
        mv "${sampleId}.cons.fa" "${sampleId}.consensus.fasta" || exit 1
        """
}

// BLAST processes (unchanged from original sangerFlow)
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
        # Ensure we have a BLAST database with proper metadata
        if [ ! -d "${db}/" ]; then
            echo "ERROR: BLAST database directory ${db} does not exist!"
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
            echo "ERROR: BLAST results file is empty or missing!"
            exit 1
        fi
        
        head -n ${params.topHits} ${query} > ${query.baseName}_topHits.tsv || exit 1
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
        if [ ${combinedFiles.size()} -eq 0 ]; then
            echo "ERROR: No files to concatenate!"
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

    // Validate samples.tsv format first (fails fast)
    validateSamples(params.samples_tsv)
    
    // Read samples.tsv: sampleId, read1 (forward SCF), read2 (reverse SCF)
    def samplesChannel = channel.fromPath( params.samples_tsv )
        | splitCsv( header: true, sep: '\t' )
        // Validate that each row has 3 columns and sampleId is not null
        | filter { it.size() == 3 && it[0] != null }
        | map { row -> tuple( row.sampleId, file(row.read1), file(row.read2) ) }

    tracyAssemble( samplesChannel )
    
    // Uncomment next lines when BLAST database is available:
    // blastn( tracyAssemble.out.consensus, params.db_path )
    // blastnTopHits(blastn.out.blastHits)
    // blastHits_concat(blastnTopHits.out.topHits.collect())
}

// Validation process to check samples.tsv format - FAILS FAST!
process validateSamples {
    tag "validation"
    
    script:
        """
        # Read the first line of the TSV file and verify headers
        header=$(head -n 1 $params.samples_tsv)
        
        # Check if the header matches exactly (tab-separated, no extra spaces)
        if [[ "$header" != "sampleId\tread1\tread2" ]]; then
            echo "ERROR: samples.tsv header must be exactly: sampleId\\tread1\\tread2"
            echo "Current header: '$header'"
            exit 1
        fi
        
        # Check if the file has data rows
        line_count=$(wc -l < $params.samples_tsv)
        if [ "$line_count" -lt 2 ]; then
            echo "ERROR: samples.tsv must have at least one data row after the header"
            exit 1
        fi
        
        # Check that all files referenced in samples.tsv exist
        tail -n +2 $params.samples_tsv | while IFS=$'\t' read -r sampleId read1 read2; do
            if [ ! -f "$read1" ]; then
                echo "ERROR: Forward file not found for sample $sampleId: $read1"
                exit 1
            fi
            if [ ! -f "$read2" ]; then
                echo "ERROR: Reverse file not found for sample $sampleId: $read2"
                exit 1
            fi
        done
        
        echo "samples.tsv format and all files validated successfully!"
        """
}

// Pre-flight check for container runtime availability - FAILS FAST!
process preFlightCheck {
    when: false  // Run this but don't emit anything, just fail if conditions not met
    
    script:
        """
        # Check Docker or Singularity availability with timeout
        echo "Checking container runtime..."
        
        if command -v docker &> /dev/null; then
            echo "✓ Docker is available"
            docker --version | head -1
        elif command -v singularity &> /dev/null; then
            echo "✓ Singularity is available"
            singularity version 2>&1 | head -1
        else
            echo "✗ ERROR: Neither Docker nor Singularity is installed!"
            echo "Please install a container runtime before running this pipeline."
            exit 1
        fi
        
        # Also check for wget/curl for Tracy download
        if command -v wget &> /dev/null || command -v curl &> /dev/null; then
            echo "✓ wget or curl is available"
        else
            echo "✗ ERROR: Neither wget nor curl is available!"
            echo "Please install one of these tools to download Tracy."
            exit 1
        fi
        
        echo "Pre-flight check passed! Pipeline can proceed."
        """
}
