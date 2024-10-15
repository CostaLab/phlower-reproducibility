library(Seurat)
library(glue)
library(dplyr)
library(glue)
library(ggplot2)
library(future.apply)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000* 1024^2)



i=1:4
samples <- glue("data/S{i}/outs/filtered_feature_bc_matrix/")
mtx_list <- lapply(samples, function(x) Read10X(x))
gene_mtx_list <- lapply(mtx_list, function(x)x[["Gene Expression"]])

gene_mtx_list <- lapply(1:length(gene_mtx_list), function(x){
                                  colnames(gene_mtx_list[[x]]) <- glue("S{x}#{colnames(gene_mtx_list[[x]])}")
                                  gene_mtx_list[[x]]})

scrna_list <- lapply(gene_mtx_list, function(x) CreateSeuratObject(counts = x,
                                                                   min.cells = 0,
                                                                   min.features=0))

scrna <- merge(x = scrna_list[[1]], y = unlist(scrna_list[2:length(scrna_list)]),
               merge.data = FALSE, project = "multiomic")
dir.create("save")
saveRDS(scrna, "save/scrna.rds")
