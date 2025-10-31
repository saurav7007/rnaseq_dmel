#!/usr/bin/env nextflow

// Include the required process
include { MERGE_COUNTS } from './processes/deg/merge_counts.nf'

// User Parameters
params.countfolder = "./results/counts"
params.outdir = "./results"

// Pipeline Logic
workflow {
    merge_counts_ch = Channel.fromPath("${params.countfolder}/*.counts.txt")
                        .collect()

    // Step 6: Merge count data into count file
    MERGE_COUNTS(merge_counts_ch)
}