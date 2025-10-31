#!/usr/bin/env nextflow

// Include the required process
include { FASTP_TRIM } from './processes/trim/fastp.nf'
include { HISAT2_INDEX; HISAT2_ALIGN } from './processes/align/hisat2.nf'
include { SAMTOOLS_SORT_INDEX } from './processes/align/samtools.nf'
include { FEATURECOUNTS } from './processes/counts/subread.nf'

// User Parameters
params.samplesheet = "./samplesheet.csv"
params.genome_fa  = "./ref/dmel.fa"
params.annotation_gtf = "./ref/dmel.gtf"
params.outdir = "./results"

// Pipeline Logic
workflow {

    // Channels
    genome_ch = Channel.fromPath(params.genome_fa)
    annot_ch = Channel.fromPath(params.annotation_gtf)
    samples_ch = Channel
                    .fromPath(params.samplesheet)
                    .splitCsv(header: true)
                    .map { row ->
                        // Always add R1 and default read_type to single-end 
                        def reads = [ file(row.read1) ]
                        def read_type = "SE"

                        // If R2 is present, add it and switch to paired-end      
                        if( row.read2?.trim() ) {
                            reads << file(row.read2)
                            read_type = "PE"
                        }

                        tuple(row.sample_id, row.condition, read_type, reads)
                    }

    // Step 1: QC & Trimming
    trimmed_ch = FASTP_TRIM(samples_ch)
                    .map { sample_id, cond, read_type, fastqs, html, json ->
                        tuple(sample_id, cond, read_type, fastqs)
                    }

    // Step 2: Build HISAT2 index
    index_ch = HISAT2_INDEX(genome_ch)
                    .collect()

    // Step 3: Align trimmed reads (.sam)
    aligned_sam_ch = HISAT2_ALIGN(trimmed_ch, index_ch)

    // Step 4: Convert .sam → sorted & indexed BAM
    final_bam_ch = SAMTOOLS_SORT_INDEX(aligned_sam_ch)
    
    // Step 5: Gene-level quantification
    counts_ch = final_bam_ch
                    .combine(annot_ch)
                    | FEATURECOUNTS
    
}
