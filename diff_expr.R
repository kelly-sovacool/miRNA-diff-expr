# Differential expression analysis with DESeq2

library(edgeR)
library(DESeq2)
library(ggplot2)
library(scales)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(ReportingTools)
library(pcaExplorer)
setwd("/Users/kelly/gdrive_umich/tewari_lab/U01_exRNA_Healthy_Persons_and_Feeding/miRNA-diff-expr")

run_deseq <- function(samples, counts, annotation_columns) {
  # Run DESeq2 pipeline
  rownames(counts) = counts$X
  rownames(samples) = samples$X
  counts$X <- NULL
  samples$X <- NULL
  samples$Participant.ID <- factor(samples$Participant.ID)
  counts <- counts[,unique(rownames(samples))]
  all(colnames(counts) == rownames(samples))
  deseq2Data <- DESeqDataSetFromMatrix(countData = counts, design = ~ Participant.ID, colData = samples)
  deseq2Data <- DESeq(deseq2Data) # this runs the whole DESeq pipeline--the individual steps can be explicitly called as well. Note that it just adds to the dds object. The counts are still kept in a specific slot in the dds object. 
  results <- results(deseq2Data, cooksCutoff=FALSE, independentFiltering=FALSE) # this extracts the results from the output
  #results_df = as.data.frame(results)
  
  # generate a report
  #report <- HTMLReport(shortName='plasma_only', title="Plasma Healthy Ctrls", basePath="/Users/kelly/gdrive_umich/tewari_lab/U01_exRNA_Healthy_Persons_and_Feeding/kelly_workdir", reportDirectory = './reports')
  #publish(deseq2Data, report, pvalueCutoff=0.1, factor = colData(deseq2Data)$Source)
  #finish(report)
  #pcaExplorer(dds=deseq2Data)
  
  # transform the data
  rlogData <- rlog(deseq2Data, blind=FALSE)
  vstData <- varianceStabilizingTransformation(deseq2Data, blind=FALSE)
  ntData <- normTransform(deseq2Data)
  
  # heatmap of transformed expression data
  select <- order(rowMeans(counts(deseq2Data,normalized=TRUE)), decreasing=TRUE)
  Source <- colData(deseq2Data)[,annotation_columns]
  annotations <- as.data.frame(Source)
  rownames(annotations) <- colnames(deseq2Data)
  pheatmap(assay(vstData)[select,], cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=annotations)
  
  # heatmap of sample-to-sample distances
  sampleDists <- dist(t(assay(rlogData)))
}


# matched plasma & serum
matched_samples <- read.csv('matched_plasma-serum_samples.csv')
matched_mir_counts <- read.csv('matched_plasma-serum_mir_counts.csv')
run_deseq(matched_samples, matched_mir_counts, c("Source"))

# plasma only
plasma_samples <- read.csv('plasma_samples.csv')
plasma_mir_counts <- read.csv('plasma_mir_counts.csv')

save.image('diff_expr_healthyCtrls.RData')
