## import pandas as pd

## ------------------------------------------------------------------------
samples <- read.csv('data/matched_plasma-serum_samples.csv')
samples <- subset(samples, Library.Generation.Set != "SetRecheck" )  # exclude SetRecheck samples
# sample X11 is outlier -- remove it and its match
outlier_id <- samples[samples$X=="X11",]$Participant.ID
samples <- samples[samples$Participant.ID != outlier_id,]
counts <- read.csv('data/matched_plasma-serum_mir_counts.csv')
samples$Participant.ID <- factor(samples$Participant.ID)  # ID is categorical, not numerical
samples$Sample.ID <- factor(samples$Sample.ID)
rownames(counts) = counts$X
rownames(samples) = samples$X
counts$X <- NULL  # remove extra column
counts <- counts[,rownames(samples)]
head(samples)
head(counts)

## ------------------------------------------------------------------------
library(edgeR)
design <- model.matrix(~Participant.ID + Source, samples)
dge = DGEList(counts = counts, samples = samples)
# require miRNAs to have CPM > 1 in at least 2 samples
countsPerMillion <- edgeR::cpm(dge)
countCheck <- countsPerMillion > 1
head(countCheck)
keep <- which(rowSums(countCheck) >= 2) 
dge <- dge[keep,]

## ------------------------------------------------------------------------
library(SingleCellExperiment)
library(scater)
reads_sce <- SingleCellExperiment(assays=list(counts=dge$counts),  colData=dge$samples)
# remove unexpressed miRNAs
keep_feature <- rowSums(counts(reads_sce) > 0) > 0
reads_sce <- reads_sce[keep_feature, ]
reads_sce <- calculateQCMetrics(reads_sce)
# create separate plasma & serum SCEs
plasma_samples <- subset(colData(reads_sce), Source == "Plasma")
plasma_counts <- as.data.frame(counts(reads_sce))[rownames(plasma_samples)]
plasma_sce <- SingleCellExperiment(assays=list(counts=as.matrix(plasma_counts)), colData=plasma_samples)
plasma_sce <- calculateQCMetrics(plasma_sce)
serum_samples <- subset(colData(reads_sce), Source == "Serum")
serum_counts <- as.data.frame(counts(reads_sce))[rownames(serum_samples)]
serum_sce <- SingleCellExperiment(assays=list(counts=as.matrix(serum_counts)), colData=serum_samples)
serum_sce <- calculateQCMetrics(serum_sce)
# log transform
cpm(reads_sce) <- calculateCPM(reads_sce)
reads_sce <- normalize(reads_sce)
logcounts(reads_sce) <- log2(calculateCPM(reads_sce) + 1)
# visualize
hist(reads_sce$total_counts, breaks=100)  # counts per sample
hist(reads_sce$total_features, breaks=100)  # counts per miRNA
plotQC(reads_sce, type = "highest-expression") + ggtitle("Plasma & Serum")
plotQC(plasma_sce, type = "highest-expression") + ggtitle("Plasma")
plotQC(serum_sce, type = "highest-expression") + ggtitle("Serum")
plotQC(reads_sce, type="explanatory-variables", variables=c("Index", "Participant.ID", "Collection.Date", "Library.Generation.Set", "MiSeq.QC.Run", "Source", 'total_features', "Age", "Race", "Sex"))

## ------------------------------------------------------------------------
plotPCA(reads_sce, exprs_values = "logcounts", colour_by = "Library.Generation.Set", size_by = "total_features", shape_by = "Source")
for (var in c("total_features", "Age", "Library.Generation.Set", "Index", "Participant.ID", "Source", "Race", "Gender", "Collection.Date")) {
  print(
    plotQC(reads_sce, type = "find-pcs", exprs_values = "logcounts", variable = var)
    )  + ggtitle(var) 
  }

## ------------------------------------------------------------------------
library(RUVSeq)
library(ggplot2)
library(mvoutlier)
dge <- calcNormFactors(dge)
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)
res <- residuals(fit, type="deviance")
set <- newSeqExpressionSet(dge$counts, phenoData=dge$samples)
ruvr_sets <- list()
for(k in 1:5) {
  ruvr_sets[[k]] <- RUVr(set, row.names(dge), k=k, res)
  assay(reads_sce, paste("RUVr k=", toString(k))) <- log2(t(t(assayData(ruvr_sets[[k]])$normalizedCounts) / colSums(assayData(ruvr_sets[[k]])$normalizedCounts) * 1e6) + 1)
}
for(n in assayNames(reads_sce)) {
  print(
        plotPCA(
            reads_sce,
            colour_by = "Library.Generation.Set",
            size_by = "total_features",
            shape_by = "Source",
            exprs_values = n
        ) + ggtitle(n)
  )
}

## ------------------------------------------------------------------------
reads_sce <- runPCA(reads_sce, use_coldata = TRUE, detect_outliers = TRUE)
outliers <- colnames(reads_sce)[reads_sce$outlier]
head(outliers)

## ------------------------------------------------------------------------
for (var in c("total_features", "Age", "Library.Generation.Set", "Index", "Participant.ID", "Source", "Race", "Gender", "Collection.Date")) {
  print(
    plotQC(reads_sce, type = "find-pcs", exprs_values = "RUVr k= 2", variable = var)
    )  + ggtitle(var) 
}

## ------------------------------------------------------------------------
library(RColorBrewer)
library(reshape2)
# Ryan's code, modified
norm.expr.matr <- exprs(reads_sce)
# Rank the mean expression values for plasma/serum miRs. Highest expression = 1
mean.expr.plasma.rank <- rank(-1*rowMeans(norm.expr.matr[, colData(reads_sce)$Source=="Plasma"]))
mean.expr.serum.rank <- rank(-1*rowMeans(norm.expr.matr[, colData(reads_sce)$Source=="Serum"]))
top_N <- 18  # top_N=18 resuls in 20 miRNAs in the plot
top.miRs <- row.names(norm.expr.matr)[mean.expr.plasma.rank <= top_N | mean.expr.serum.rank <= top_N] # get the names of the top miRs in plasma or serum
norm.expr.top <- norm.expr.matr[top.miRs, ] # Get the expression matrix for the top mIRs
norm.expr.melt <- reshape2::melt(norm.expr.top) # Convert your normalized expression matrix to a 3 column data.frame (row.name, col.name, expression value). Melt is in the dpylr package, I believe.
colnames(norm.expr.melt) <- c("miR.ID", "MT.Unique.ID", "norm.expr") # just so it's easier for me to tell you which columns I'm using.
norm.expr.melt$Source <- ""
for (row_num in 1:nrow(norm.expr.melt)){  # Pull the plasma/serum source value from the column metadata.
  mt_unique_id <- norm.expr.melt[row_num, ]$MT.Unique.ID
  norm.expr.melt[row_num, "Source"] <- as.character(samples[mt_unique_id,"Source"])
}
# This would be a simple boxplot with plasma/serum side-by-side. Overlaying the boxes over points takes a little more tweaking to get the dodge/width right, but it's doable.
ggplot(norm.expr.melt, aes(x=reorder(miR.ID, norm.expr, FUN=median), y=norm.expr, fill=Source)) + geom_boxplot(pos="dodge", outlier.size=0.5) + 
  ggtitle("Top 20 Expressed miRNAs") + ylab("Normalized Expression") + xlab("miRNA ID") + 
  theme(panel.grid.major.x = element_line(size = 0.5, linetype = 'solid', colour = "grey"), 
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(norm.expr.melt, aes(x=reorder(miR.ID, norm.expr, FUN=median), y=norm.expr, fill=Source)) + geom_boxplot(pos="dodge", outlier.size=0.5) + 
  ggtitle("Top 20 Expressed miRNAs") + ylab("Normalized Expression") + xlab("miRNA ID") + 
  theme(panel.grid.major.y = element_line(size = 0.5, linetype = 'solid', colour = "grey"), 
        panel.grid.major.x = element_blank()) + coord_flip()

## ------------------------------------------------------------------------
library(pheatmap)
# based on Ryan's code, pass RUVSeq-corrected data to edgeR
ruvr2 <- ruvr_sets[[2]]
design <- model.matrix(~ Participant.ID + Source + W_1 + W_2, pData(ruvr2))
dge <- estimateDisp(dge, design = design, tagwise = TRUE, robust = TRUE)
fit <- glmFit(dge, design)
# contrast source (plasma vs serum)
lrt <- glmLRT(fit, coef="SourceSerum")
results <- data.frame(topTags(lrt, n=Inf, sort.by="PValue", p.value=0.05))
sig_miRs = list()
for (logFC_threshold in 1:2) {
  sig_miRs[[logFC_threshold]] <- rownames(results[results$PValue < 0.05 & abs(results$logFC) > logFC_threshold,])  # filter by p-value and logFC
}
sig_miRs[[3]] <- rownames(results[results$PValue < 0.05,][0:50,])
# heatmaps of significant DE miRNAs
Source <- pData(ruvr2)[,c("Source")]
annotations <- as.data.frame(Source)
rownames(annotations) <- rownames(pData(ruvr2))
for(miR_list in sig_miRs) {
  pheatmap(norm.expr.matr[miR_list,], cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=annotations, main="Differentially Expressed miRNAs")
}

## ------------------------------------------------------------------------
library(DESeq2)  # TODO: use filtered & RUVSeq-corrected data
deseq2Data <- DESeqDataSetFromMatrix(countData = dge$counts, design = ~ Participant.ID + Source, colData = dge$samples)
deseq2Data <- DESeq(deseq2Data)
results <- results(deseq2Data, cooksCutoff=FALSE, independentFiltering=FALSE)
head(as.data.frame(results))

## ------------------------------------------------------------------------
   # TODO: only use top DE miRNAs
vstData <- varianceStabilizingTransformation(deseq2Data, blind=FALSE)  # transform the data
select <- order(rowMeans(counts(deseq2Data,normalized=TRUE)), decreasing=TRUE)
Source <- colData(deseq2Data)[,c("Source")]
annotations <- as.data.frame(Source)
rownames(annotations) <- colnames(deseq2Data)
pheatmap(assay(vstData)[select,], cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=annotations)

## ------------------------------------------------------------------------
sampleDists <- dist(t(assay(vstData)))  # compute sample-to-sample distances
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vstData$Participant.ID, vstData$Source, sep=" - ")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

## ------------------------------------------------------------------------
plotPCA(
            reads_sce,
            colour_by = "Age",
            shape_by = "Gender",
            size_by = "total_features",
            exprs_values = "RUVr k= 2"
) + ggtitle("RUVSeq-Normalized Expression (k=2)") 

## ------------------------------------------------------------------------
save.image("diff_expr_plasma_vs_serum.RData")

