process SAMTOOLS_SORT_INDEX {
    tag { sample_id }
    publishDir "${params.outdir}/aligned", mode: 'copy'
    container "quay.io/biocontainers/samtools:1.22.1--h96c455f_0"

    input:
    tuple val(sample_id),
          val(cond),
          path(sam_file)

    output:
    tuple val(sample_id),
          val(cond),
          path("${sample_id}.sorted.bam")

    script:
    """
    samtools view -bS ${sam_file} | samtools sort -@ 8 -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}
