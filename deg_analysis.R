# Install necessary packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("DEGreport")
install.packages("pheatmap")
install.packages("tidyverse")

# Calling necessary libraries
library("tidyverse")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("apeglm")
library("DEGreport")

# Creating sample information (Metadata)
samples_info <- read.csv("samplesheet.csv", row.names = 1) %>%
  mutate(
    #sample_id = paste0(sample_id, ".sorted.bam"),
    record_type = if_else(nzchar(trimws(read2)), "pair-end", "single-end")
  ) %>%
  select(condition, batch, record_type) %>%
  { `rownames<-`(., paste0(rownames(.), ".sorted.bam")) }

# Reading the merged count data (Count data)
count_data <- read.csv(
                file="results/counts/all.counts.txt",
                sep = "\t",
                header = TRUE,
                row.names = 1
              )

# Reorder Count columns based on row order of metadata (sample_info)
count_data <- count_data[, rownames(samples_info)]

# Construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = samples_info,
                              design = ~ batch + condition
                              )

# Keep rows with atleast 10 reads across samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Specifying the reference level
dds$condition <- relevel(dds$condition, ref = "WT")

# Run DESeq
dds <- DESeq(dds)

# To check if the sames are clustered based on bactch
vsd <- vst(dds)
plotPCA(vsd, intgroup = c("batch", "condition"))

# Results
res <- results(dds)
resultsNames(dds)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20] 

df <- as.data.frame(colData(dds)[,c("condition","record_type")])

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$record_type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup = c("condition", "record_type"))