# Ryan's script for getting Canonical miRNA counts
library(data.table)
library(Hmisc)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(stringr)
setwd("/Users/kelly/gdrive_umich/tewari_lab/U01_exRNA_Healthy_Persons_and_Feeding/kelly_workdir")

# Import unfiltered counts. We will do our own filtering
counts.files <- Cs(isoform.Counts.csv, isoform.Counts.Mismatch.csv,  miR.Counts.csv, miR.Counts.Mismatch.csv)

# Import sample sheets and fix names/formatting
sample.annot <- fread("sample_sheet.csv")
sample.annot <- subset(sample.annot, MISEQ.QC.PASS=="PASS")
sample.annot[Library.Generation.Set=="SetRecheck", MiSeq.QC.Run:="Run9"]

sample.annot$Gender <- tolower(sample.annot$Gender)
sample.annot$Race <- gsub("/", "", gsub(" ", ".", tolower(sample.annot$Race)))
sample.annot$Library.Generation.Set <- gsub(" ", "", sample.annot$Library.Generation.Set)
sample.annot$Sample.ID <- gsub("\\*$", "", sample.annot$Sample.ID)
sample.annot$Participant.ID <- gsub("\\*$", "", sample.annot$Participant.ID)
sample.annot[, MT.Unique.ID:=factor(paste0("X", MT.Unique.ID), levels=paste0("X", 1:max(MT.Unique.ID)))]
setkey(sample.annot, MT.Unique.ID)

group_directories <- "/Users/kelly/gdrive_umich/tewari_lab/U01 exRNA Healthy Persons and Feeding/processed_data/201704_U01_HiSeq_HealthyControlStudy/02-map_contaminates/output"

# Function imports files and assigns them to variables with the name of the file minus the ".csv" suffix
import.files <- function(in.file, group_dirs=group_directories, sample_annotations=sample.annot){
  in.files <- dir(group_dirs, pattern=in.file, recursive=TRUE, full.names=TRUE) # searches current directory recursively under group_dirs for pattern in.file and returns path to files
  var.name <- tolower(sub(".csv", "", in.file)) # drops the .csv suffix to generate names for variables that will be assigned in.GlobalEnv. Willl export as DGEList, so adds dge
  var.name.long <- paste0(var.name, ".long.dt") # creates a second set of variable names for the data table left in melted form 
  var.name <- paste0(var.name, ".dge")
  message(paste0("# Reading ", in.file, " from group directories")) # prints message
  # Iterates over the files in "in.files". 
  # For each file, imports as a data.table, then melts to a data table. Melted data tables have the first two columns as the original row and column ids, respectively. 
  # The third column is the counts originally in the wide matrix
  # The melted tables are concatenated across the input files
  dt.long <- do.call(rbind, 
                     sapply(in.files, 
                            USE.NAMES=FALSE,
                            simplify=FALSE,
                            FUN=function(x){
                              dt <- fread(x, header=TRUE)
                              dt.m <- melt.data.table(dt, id.vars=1)
                            }))
  dt.cast <- dcast.data.table(dt.long, ... ~ variable, value.var="value", fun.aggregate=sum, fill=0) # re-casts the concatenated data table, filling missing values with 0
  dt.cast.matr <- as.matrix(data.frame(dt.cast, row.names=1)) # converts casted data table to matrix with row and column names
  colnames(dt.cast.matr) <- make.names(strsplit2(colnames(dt.cast.matr), split="_")[,1]) # Converts the column names to drop the full file name and leave just the MT.Unique.ID
  
  setkey(dt.long, variable)
  id.dt <- dt.long[, .N, by=variable][, .(MT.Unique.ID=paste0("X", tstrsplit(variable, split="_")[1])), by=variable]
  dt.long2 <- merge(dt.long, id.dt)
  dt.long2[, variable:=NULL]
  setnames(dt.long2, "value", "count")
  dt.long2 <- subset(dt.long2, count>0)
  dt.long2[, MT.Unique.ID:=factor(MT.Unique.ID, levels=levels(sample_annotations$MT.Unique.ID))]
  if(str_detect(var.name, "isoform")){
    mir.id.split <- data.table(cbind(row.names(dt.cast.matr), strsplit2(row.names(dt.cast.matr), split="_")))
    setnames(mir.id.split, c("miR.Isoform.ID", "miR.ID", "Offset.5p.tmp", "Offset.3p.tmp"))
    mir.id.split[, `:=`(Offset.5p=as.numeric(Offset.5p.tmp), Offset.3p=as.numeric(Offset.3p.tmp))][, c("Offset.5p.tmp", "Offset.3p.tmp"):=NULL]
    setkey(mir.id.split, miR.Isoform.ID)
    setnames(dt.long2, 1, "miR.Isoform.ID")
    setkey(dt.long2, miR.Isoform.ID)
    dt.long3 <- merge(dt.long2, mir.id.split)
    setkey(dt.long3, MT.Unique.ID)
    setkey(sample_annotations, MT.Unique.ID)
    dt.long3 <- merge(sample_annotations, dt.long3)
    assign(var.name.long, dt.long3, envir=.GlobalEnv)
    mir.id.split.df <- data.frame(mir.id.split, row.names=1)
    sample.annot.df <- data.frame(sample_annotations, row.names=1)[colnames(dt.cast.matr),]
    dt.cast.dge <- DGEList(dt.cast.matr, genes=mir.id.split.df)
    dt.cast.dge$samples <- cbind(dt.cast.dge$samples, sample.annot.df)
    assign(var.name, dt.cast.dge, envir=.GlobalEnv)  # assigns the resulting count matrix to the variable derived from the file name
  } else{
    setkey(dt.long2, MT.Unique.ID)
    setkey(sample_annotations, MT.Unique.ID)
    dt.long2 <- merge(sample_annotations, dt.long2)
    dt.cast.dge <- DGEList(dt.cast.matr)
    sample.annot.df <- data.frame(sample_annotations, row.names=1)[colnames(dt.cast.matr),]
    dt.cast.dge$samples <- cbind(dt.cast.dge$samples, sample.annot.df)
    assign(var.name, dt.cast.dge, envir=.GlobalEnv)
    assign(var.name.long, dt.long2, envir=.GlobalEnv)
  } 
  message(paste0("# Done with ", in.file, ". Output Written assigned to ", var.name, "(DGEList) and ", var.name.long, " (melted data table).")) # prints progress message
}

sapply(counts.files, import.files) # Iterates import.files function over counts.files
# Reading isoform.Counts.csv from group directories
# Done with isoform.Counts.csv. Output Written assigned to isoform.Counts.dge(DGEList) and isoform.counts.long.dt (melted data table).
# Reading isoform.Counts.Mismatch.csv from group directories
# Done with isoform.Counts.Mismatch.csv. Output Written assigned to isoform.Counts.Mismatch.dge(DGEList) and isoform.counts.mismatch.long.dt (melted data table).
# Reading miR.Counts.csv from group directories
# Done with miR.Counts.csv. Output Written assigned to miR.Counts.dge(DGEList) and mir.counts.long.dt (melted data table).
# Reading miR.Counts.Mismatch.csv from group directories
# Done with miR.Counts.Mismatch.csv. Output Written assigned to miR.Counts.Mismatch.dge(DGEList) and mir.counts.mismatch.long.dt (melted data table).

# do some isoform counts filtering
studies <- "Healthy Controls"
setkey(isoform.counts.long.dt, miR.ID)
isoform.counts.long.dt.feeding.study <- subset(isoform.counts.long.dt, Study%in%studies)
isoform.counts.long.dt.feeding.study[, `:=`(abs.offset.5p=abs(Offset.5p), abs.offset.3p=abs(Offset.3p), tot.abs.offset=abs(Offset.5p)+abs(Offset.3p))]
isoform.counts.long.dt.feeding.study[, total.counts.sample.unfiltered:=sum(count), by=MT.Unique.ID]
max.offset.5p=1
max.offset.3p=2
max.total.offset=2
by.cols <- c(colnames(sample.annot), c("miR.ID", "total.counts.sample.unfiltered"))
canonical.isoform.counts.long.dt.feeding.study <- isoform.counts.long.dt.feeding.study[, .(tot.count.miR.sample=sum(count), tot.count.canonical.miR.sample=sum(ifelse(abs.offset.5p<=max.offset.5p & abs.offset.3p<=max.offset.3p & tot.abs.offset<=max.total.offset, count, 0L))),  by=by.cols]
canonical.isoform.counts.long.dt.feeding.study[, `:=`(cpm.all.unfiltered=(tot.count.miR.sample*10^6)/total.counts.sample.unfiltered,
                                                      cpm.canonical.unfiltered=(tot.count.canonical.miR.sample*10^6)/total.counts.sample.unfiltered ,
                                                      percent.canonical=tot.count.canonical.miR.sample/tot.count.miR.sample)]
canonical.isoform.counts.long.dt.feeding.study[, n.samples.with.canonical.expr:=sum(ifelse(tot.count.canonical.miR.sample>0, 1, 0)), by=miR.ID]
total.samples.feeding.study <- dim(canonical.isoform.counts.long.dt.feeding.study[, .N, by=MT.Unique.ID])[1]

canonical.isoform.counts.long.dt.feeding.study.summary <- canonical.isoform.counts.long.dt.feeding.study[, .(n.samples.with.expr=sum(ifelse(tot.count.miR.sample>0, 1, 0)),
                                                                                                             n.samples.with.canonical.expr=sum(ifelse(tot.count.canonical.miR.sample>0, 1, 0)), 
                                                                                                             mean.percent.canonical=mean(percent.canonical),
                                                                                                             sum.cpm.canonical=sum(cpm.canonical.unfiltered),
                                                                                                             sum.cpm.all=sum(cpm.all.unfiltered)), by=miR.ID]
canonical.isoform.counts.long.dt.feeding.study.summary[, `:=`(mean.cpm.canonical=ifelse(n.samples.with.canonical.expr==0, NA, sum.cpm.canonical/n.samples.with.canonical.expr),
                                                              mean.cpm.all=sum.cpm.all/n.samples.with.expr)]

#
canonical.isoform.counts.long.dt.feeding.study.cast <- dcast.data.table(canonical.isoform.counts.long.dt.feeding.study, miR.ID ~ MT.Unique.ID, value.var="tot.count.canonical.miR.sample", fill=0, fun.aggregate=sum)
canonical.counts.miR.dge <- DGEList(data.frame(canonical.isoform.counts.long.dt.feeding.study.cast, row.names=1))
canonical.counts.miR.dge$samples <- merge(canonical.counts.miR.dge$samples, mir.counts.dge$samples[colnames(canonical.counts.miR.dge),], by=0)
canonical.counts.miR.dge$samples$timepoint <- strsplit2(canonical.counts.miR.dge$samples$Sample.ID, split="-")[,2]
colnames(canonical.counts.miR.dge$samples) <- sub(".x$", "", colnames(canonical.counts.miR.dge$samples))

save.image("Canonical_Counts_miRs_HiSeq_Healthy_Controls.RData")

