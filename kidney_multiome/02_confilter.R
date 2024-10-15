library(Seurat)
library(glue)
library(dplyr)
library(glue)
library(ggplot2)
library(future.apply)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000* 1024^2)

meta <- read.csv("save/cellColData.tsv", sep="\t")
scrna <- readRDS("save/scrna.rds")

nms <- intersect(rownames(meta), colnames(scrna))
message("#cells: ", length(nms))
scrna <- subset(scrna, cells=nms)
meta <- meta[nms,]

scrna@meta.data <- cbind(scrna@meta.data, meta)
scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^mt-|^MT-")
scrna[["percent.ribo"]] <- PercentageFeatureSet(scrna, pattern = "^Rpl|^Rps|^RPL|^RPS")
scrna <- subset(scrna, subset = nFeature_RNA > 400 &  nCount_RNA < 400000 & percent.mt < 5)
message("dim: ", dim(scrna))

## save confiltered cells metadata
write.table(scrna@meta.data,
                file="save/cellColData_confiltered.tsv",
                sep="\t",
                quote=FALSE,
                row.names=TRUE)

saveRDS(scrna, "save/scrna_confiltered.rds")


