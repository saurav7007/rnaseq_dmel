/*
 * RNA-Seq pipeline: Step 6: Merge count data into count file
 * Input : All the count files from FEATURECOUNTS
 * Output: Single merged count file
 */

process MERGE_COUNTS {
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
        path(count_files)

    output:
        path("all.counts.txt")

    script:
    '''
    # Extract the counts per sample
    for file in *.counts.txt; do
        name=$(basename --suffix .counts.txt ${file})
        sed "1d" ${file} | cut -f7 > ${name}_counts_clean.txt;
    done

    # Extract gene_id
    ls -1 *.counts.txt | head -1 | xargs cut -f1 | tail -n +2 > genes.txt

    # Merge all the count data into one file
    paste genes.txt *counts_clean.txt > all.counts.txt
    '''
}