/*
 * RNA-Seq pipeline — Step 2: Alignment with HISAT2
 * Input: trimmed FASTQ pairs  (from fastp)
 * Output: sorted BAM
 */

/*
 * HISAT2_INDEX — build genome index once
 * Input: reference genome FASTA
 * Output: HISAT2 index files (*.ht2)
 */

process HISAT2_INDEX {
    publishDir "${params.outdir}/index", mode: 'copy'
    container 'quay.io/biocontainers/hisat2:2.2.1--h503566f_8'

    input:
    path genome_fa

    output:
    path "dmel_index.*.ht2"

    script:
    """
    echo "Building HISAT2 index from ${genome_fa}..."
    hisat2-build --threads 8 ${genome_fa} dmel_index
    """
}

process HISAT2_ALIGN {
    tag { sample_id }
    //publishDir "${params.outdir}/aligned", mode: 'copy'
    container 'quay.io/biocontainers/hisat2:2.2.1--h503566f_8'

    input:
    tuple val(sample_id),
          val(cond),
          val(read_type),
          path(reads)
    path index_files

    output:
    tuple val(sample_id),
          val(cond),
          val(read_type),
          path("${sample_id}.sam")

    script:
    """
    if [[ "${read_type}" == "PE" ]]; then
        # Paired-end
        hisat2 -x dmel_index \
            -1 ${reads[0]} -2 ${reads[1]} \
            -S ${sample_id}.sam \
            --threads 4
    else
        # Single-end
        hisat2 -x dmel_index \
           -U ${reads[0]} \
           -S ${sample_id}.sam \
           --threads 4
    fi
    """
}
