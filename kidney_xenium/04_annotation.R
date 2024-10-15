library(Seurat)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(dplyr)
library(plyr)
library(ggpubr)
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggpp))
suppressPackageStartupMessages(library(ggalluvial))
mute = function(x) suppressWarnings(suppressMessages(x))


source("util/util.R")


object <- readRDS("save/all_353_pca_cleaned.Rds")
Idents(object) <- "seurat_clusters"

annotation = annotation_353

annotation <- plyr::mapvalues(as.character(Idents(object)), from = names(annotation), to = annotation)
object$annotation <- annotation

write.csv(object@meta.data, "save/all_353_cleaned_meta.csv")
saveRDS(object, "save/all_353_cleaned_pca_annotation.Rds")
