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

ids =c("353_df",
       "353_si",
       "353_d12",
       "353_d18"
      )




if(F){
  object <- readRDS("save/all_353_cleaned_pca_annotation.Rds")
  Idents(object) <- "seurat_clusters"
  cluster_to_remove <- c("4", "10", "12", "13", "14")
  object <- subset(object, idents = cluster_to_remove, invert = TRUE)
  #remove null factor levels of seurat_clusters
  object$seurat_clusters <- droplevels(object$seurat_clusters)

  write.csv(object@meta.data, "save/07_all_353_cleaned_subset_meta.csv")
  saveRDS(object, "save/07_all_353_cleaned_pca_annotation_subset.Rds")

  pca <- Embeddings(object[["pca"]])
  write.csv(pca, "save/07_all_353_cleaned_subset_pca.csv")

  umap <- Embeddings(object[["umap"]])
  write.csv(umap, "save/07_all_353_cleaned_subset_umap.csv")




  New<- JoinLayers(object, overwrite = TRUE, assay = "Xenium")
  rna_counts <- GetAssayData(New, assay = "Xenium", slot = "counts")


  DropletUtils::write10xCounts(path="save/07_all_353_cleaned_pca_annotation_subset.h5", x=rna_counts, type="HDF5")



  ##subset each sample
  for(x in ids){## filtering by manually selected region.
    message(date(), " ", x)
    obj <- readRDS(glue::glue("save/iPS{x}_obj_cleaned.Rds"))
    cell_ids <- colnames(object)[object$sample == paste0("iPS",x)]
    ##remove cell_ids predix paste0("iPS",x)
    cell_ids <- gsub(paste0("^iPS",x, "_"), "", cell_ids)

    assertthat::assert_that(all(cell_ids %in% colnames(obj)))
    obj <- subset(obj, cells = cell_ids)
    saveRDS(obj, glue::glue("save/07_iPS{x}_obj_cleaned_subset.Rds"))
  }
}

