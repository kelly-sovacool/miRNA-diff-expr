## import pandas as pd

## ------------------------------------------------------------------------
samples <- read.csv('data/plasma_samples.csv')
counts <- read.csv('data/plasma_mir_counts.csv')
samples <- subset(samples, Library.Generation.Set != "SetRecheck" )  # exclude SetRecheck samples
samples <- samples[samples$X != "X11" & samples$X != "X92",] # remove outliers
samples$Participant.ID <- factor(samples$Participant.ID)  # ID is categorical, not numerical
rownames(counts) = counts$X
rownames(samples) = samples$X
counts$X <- NULL  # remove extra column
samples$X <- NULL
counts <- counts[,unique(rownames(samples))]
head(samples)
head(counts)

## ------------------------------------------------------------------------
library(edgeR)
design <- model.matrix(~0+Race+Gender+Age, samples)
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
# log transform
cpm(reads_sce) <- calculateCPM(reads_sce)
reads_sce <- normalize(reads_sce)
logcounts(reads_sce) <- log2(calculateCPM(reads_sce) + 1)
# visualize
hist(reads_sce$total_counts, breaks=100)  # counts per sample
hist(reads_sce$total_features, breaks=100)  # counts per miRNA
plotQC(reads_sce, type = "highest-expression")
plotQC(reads_sce, type="explanatory-variables", variables=c("Index", "Participant.ID", "Collection.Date", "Library.Generation.Set", "MiSeq.QC.Run", 'total_features', "Age", "Race", "Gender"))

## ------------------------------------------------------------------------
plotPCA(reads_sce, exprs_values = "logcounts", colour_by = "Library.Generation.Set", size_by = "total_features")
plotPCA(reads_sce, exprs_values = "logcounts", colour_by = "Race", shape_by="Gender", size_by = "total_features")
for (var in c("total_features", "Age", "Library.Generation.Set", "Index", "Participant.ID", "Race", "Gender", "Collection.Date")) {
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
            exprs_values = n
        ) + ggtitle(paste("Batch",n))
  )
  print(
        plotPCA(
            reads_sce,
            colour_by = "Race",
            shape_by = "Gender",
            size_by = "total_features",
            exprs_values = n
        ) + ggtitle(paste("Demographics", n))        
  )
}

## ------------------------------------------------------------------------
reads_sce <- runPCA(reads_sce, use_coldata = TRUE, detect_outliers = TRUE)
outliers <- colnames(reads_sce)[reads_sce$outlier]
head(outliers)

## ------------------------------------------------------------------------
for (var in c("total_features", "Age", "Library.Generation.Set", "Index", "Participant.ID", "Race", "Gender", "Collection.Date")) {
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
mean.expr.male.rank <- rank(-1*rowMeans(norm.expr.matr[, colData(reads_sce)$Gender=="Male"]))
mean.expr.female.rank <- rank(-1*rowMeans(norm.expr.matr[, colData(reads_sce)$Gender=="Female"]))
top_N <- 20
top.miRs <- row.names(norm.expr.matr)[mean.expr.male.rank <= top_N | mean.expr.female.rank <= top_N] # get the names of the top miRs in plasma or serum
norm.expr.top <- norm.expr.matr[top.miRs, ] # Get the expression matrix for the top mIRs
norm.expr.melt <- reshape2::melt(norm.expr.top) # Convert your normalized expression matrix to a 3 column data.frame (row.name, col.name, expression value). Melt is in the dpylr package, I believe.
colnames(norm.expr.melt) <- c("miR.ID", "MT.Unique.ID", "norm.expr") # just so it's easier for me to tell you which columns I'm using.
norm.expr.melt$Gender <- ""
for (row_num in 1:nrow(norm.expr.melt)){  # Pull the Gender values from the column metadata.
  mt_unique_id <- norm.expr.melt[row_num, ]$MT.Unique.ID
  norm.expr.melt[row_num, "Gender"] <- as.character(samples[mt_unique_id,"Gender"])
}
# This would be a simple boxplot with plasma/serum side-by-side. Overlaying the boxes over points takes a little more tweaking to get the dodge/width right, but it's doable.
ggplot(norm.expr.melt, aes(x=reorder(miR.ID, norm.expr, FUN=median), y=norm.expr, fill=Gender)) + geom_boxplot(pos="dodge", outlier.size=0.5) + 
  ggtitle("Top 20 Expressed miRNAs") + ylab("Normalized Expression") + xlab("miRNA ID") + 
  theme(panel.grid.major.x = element_line(size = 0.5, linetype = 'solid', colour = "grey"), 
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(norm.expr.melt, aes(x=reorder(miR.ID, norm.expr, FUN=median), y=norm.expr, fill=Gender)) + geom_boxplot(pos="dodge", outlier.size=0.5) + 
  ggtitle("Top 20 Expressed miRNAs") + ylab("Normalized Expression") + xlab("miRNA ID") + 
  theme(panel.grid.major.y = element_line(size = 0.5, linetype = 'solid', colour = "grey"), 
        panel.grid.major.x = element_blank()) + coord_flip()

## ------------------------------------------------------------------------
# Rank the mean expression values for plasma/serum miRs. Highest expression = 1
mean.expr.afr.rank <- rank(-1*rowMeans(norm.expr.matr[, colData(reads_sce)$Race=="African_American"]))
mean.expr.asn.rank <- rank(-1*rowMeans(norm.expr.matr[, colData(reads_sce)$Race=="Asian"]))
mean.expr.wht.rank <- rank(-1*rowMeans(norm.expr.matr[, colData(reads_sce)$Race=="White"]))
top_N <- 20
top.miRs <- row.names(norm.expr.matr)[mean.expr.afr.rank <= top_N | mean.expr.asn.rank <= top_N | mean.expr.wht.rank <= top_N] # get the names of the top miRs in plasma or serum
norm.expr.top <- norm.expr.matr[top.miRs, ] # Get the expression matrix for the top mIRs
norm.expr.melt <- reshape2::melt(norm.expr.top) # Convert your normalized expression matrix to a 3 column data.frame (row.name, col.name, expression value). Melt is in the dpylr package, I believe.
colnames(norm.expr.melt) <- c("miR.ID", "MT.Unique.ID", "norm.expr") # just so it's easier for me to tell you which columns I'm using.
norm.expr.melt$Race <- ""
for (row_num in 1:nrow(norm.expr.melt)){  # Pull the Gender values from the column metadata.
  mt_unique_id <- norm.expr.melt[row_num, ]$MT.Unique.ID
  norm.expr.melt[row_num, "Race"] <- as.character(samples[mt_unique_id,"Race"])
}
norm.expr.melt <- norm.expr.melt[norm.expr.melt$Race %in% c("African_American", "Asian", "White"),]
# This would be a simple boxplot with plasma/serum side-by-side. Overlaying the boxes over points takes a little more tweaking to get the dodge/width right, but it's doable.
ggplot(norm.expr.melt, aes(x=reorder(miR.ID, norm.expr, FUN=median), y=norm.expr, fill=Race)) + geom_boxplot(pos="dodge", outlier.size=0.5) + 
  ggtitle("Top 20 Expressed miRNAs") + ylab("Normalized Expression") + xlab("miRNA ID") + 
  theme(panel.grid.major.x = element_line(size = 0.5, linetype = 'solid', colour = "grey"), 
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(norm.expr.melt, aes(x=reorder(miR.ID, norm.expr, FUN=median), y=norm.expr, fill=Race)) + geom_boxplot(pos="dodge", outlier.size=0.5) + 
  ggtitle("Top 20 Expressed miRNAs") + ylab("Normalized Expression") + xlab("miRNA ID") + 
  theme(panel.grid.major.y = element_line(size = 0.5, linetype = 'solid', colour = "grey"), 
        panel.grid.major.x = element_blank()) + coord_flip()

## ------------------------------------------------------------------------
library(pheatmap)
ruvr2 <- ruvr_sets[[2]]
norm.expr.matr <- norm.expr.matr[order(rowMeans(norm.expr.matr), decreasing=TRUE),]
Demographics <- pData(ruvr2)[,c("Age","Race", "Gender")]
annotations <- as.data.frame(Demographics)
rownames(annotations) <- rownames(pData(ruvr2))
pheatmap(norm.expr.matr[0:50,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=annotations, main="Top 50 Expressed miRNAs")

## ------------------------------------------------------------------------
# based on Ryan's code, pass RUVSeq-corrected data to edgeR
ruvr2 <- ruvr_sets[[2]]
design <- model.matrix(~ 0 + Race + Gender + Age + W_1 + W_2, pData(ruvr2))
dge <- estimateDisp(dge, design = design, tagwise = TRUE, robust = TRUE)
fit <- glmFit(dge, design)

## ------------------------------------------------------------------------
lrt_gender <- glmLRT(fit, coef = "GenderMale")
results_gender <- data.frame(topTags(lrt_gender, n=Inf, sort.by="PValue", p.value=0.05))
sig_miRs_gender = list()
logFC_threshold <- 1
sig_miRs_gender[[1]] <- rownames(results_gender[results_gender$PValue < 0.05 & abs(results_gender$logFC) > logFC_threshold,])  # filter by p-value and logFC
sig_miRs_gender[[2]] <- rownames(results_gender[results_gender$PValue < 0.05,])
# heatmaps of significant DE miRNAs
Gender <- pData(ruvr2)[,c("Gender")]
annotations <- as.data.frame(Gender)
rownames(annotations) <- rownames(pData(ruvr2))
for(miR_list in sig_miRs_gender) {
  pheatmap(norm.expr.matr[miR_list,], cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=annotations, main="Differentially Expressed miRNAs")
}

## ------------------------------------------------------------------------
lrt_race <- glmLRT(fit, contrast=makeContrasts(asn=(RaceAsian-(RaceAfrican_American+RaceWhite)/2),
                                          afr=(RaceAfrican_American-(RaceAsian+RaceWhite)/2),
                                          wht=(RaceWhite-(RaceAfrican_American+RaceAsian)/2), levels=design))
results_race <- data.frame(topTags(lrt_race, n=Inf, sort.by="PValue", p.value=0.05))
sig_miRs_race = list()
sig_miRs_race[[1]] <- rownames(results_race[results_race$PValue < 0.05 & abs(results_race$logFC.asn) > logFC_threshold | abs(results_race$logFC.afr) > logFC_threshold | abs(results_race$logFC.wht) > logFC_threshold,])  # filter by p-value and logFC
sig_miRs_race[[2]] <- rownames(results_race[results_race$PValue < 0.05,])
# heatmaps of significant DE miRNAs
Demographics <- pData(ruvr2)[,c("Race", "Gender")]
annotations <- as.data.frame(Demographics)
rownames(annotations) <- rownames(pData(ruvr2))
for(miR_list in sig_miRs_race) {
  pheatmap(norm.expr.matr[miR_list,], cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=annotations, main="Differentially Expressed miRNAs")
}

## ------------------------------------------------------------------------
plotPCA(
            reads_sce,
            colour_by = "Age",
            shape_by = "Gender",
            size_by = "total_features",
            exprs_values = "RUVr k= 2"
) + ggtitle("RUVSeq-Normalized Expression (k=2)") 

## ------------------------------------------------------------------------
save.image("diff_expr_plasma_only.RData")

