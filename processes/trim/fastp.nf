/*
 * RNA-Seq pipeline — Step 1: QC and trimming with fastp
 * Input : paired-end FASTQ files
 * Output: trimmed FASTQs + HTML + JSON reports
*/

process FASTP_TRIM {
   tag { sample_id }
   publishDir "${params.outdir}/trimmed", mode: "copy"
   container "quay.io/biocontainers/fastp:1.0.1--heae3180_0"

   input:
   tuple val(sample_id),
         val(cond),
         path(reads)

   output:
   tuple val(sample_id),
         val(cond),
         path("*.trim.fastq.gz"),
         path("*.fastp.html"),
         path("*.fastp.json")

    script:
    """
    if [ ${reads.size()} -eq 2 ]; then
        # Paired-end
        fastp \
            -i ${reads[0]} \
            -I ${reads[1]} \
            -o ${sample_id}_R1.trim.fastq.gz \
            -O ${sample_id}_R2.trim.fastq.gz \
            -h ${sample_id}.fastp.html \
            -j ${sample_id}.fastp.json \
            --detect_adapter_for_pe \
            --thread 8
    else
        # Single-end
        fastp \
            -i ${reads[0]} \
            -o ${sample_id}_R1.trim.fastq.gz \
            -h ${sample_id}.fastp.html \
            -j ${sample_id}.fastp.json \
            --thread 8
    fi
    """
}
