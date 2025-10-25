/*
 * RNA-Seq pipeline — Step 5: Gene-level quantification with featureCounts
 * Input : sorted BAM files from SAMTOOLS_SORT_INDEX
 * Output: raw read-count tables
 */

process FEATURECOUNTS {
    tag { sample_id }
    publishDir "${params.outdir}/counts", mode: 'copy'
    container "quay.io/biocontainers/subread:2.1.1--h577a1d6_0"

    input:
    tuple val(sample_id),
          val(cond),
          val(read_type),
          path(bam),
          path(annot_gtf)

    output:
    tuple val(sample_id),
          val(cond),
          val(read_type),
          path("${sample_id}.counts.txt")

    script:
    """
    if [[ "${read_type}" == "PE" ]]; then
        # Paired-end
        featureCounts \
            -p \
            -T 8 \
            -a ${annot_gtf} \
            -o ${sample_id}.counts.txt \
            ${bam}
    else
        # Single-end
        featureCounts \
            -T 8 \
            -a ${annot_gtf} \
            -o ${sample_id}.counts.txt \
            ${bam}
    fi
    """
}
