library(Seurat)
library(glue)
library(dplyr)
library(glue)
library(ggplot2)
library(future.apply)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000* 1024^2)


scrna <- readRDS("save/scrna_confiltered.rds")

all.genes <- rownames(x = scrna)

s.genes <- paste0("^", cc.genes$s.genes, "$", collapse = "|")
s.genes <- all.genes[grepl(s.genes, all.genes, ignore.case = TRUE)]

g2m.genes <- paste0("^", cc.genes$g2m.genes, "$", collapse = "|")
g2m.genes <- all.genes[grepl(g2m.genes, all.genes, ignore.case = TRUE)]

scrna <- ScaleData(scrna)
scrna <- RunPCA(scrna, features = c(s.genes, g2m.genes), reduction.name="BCELLCYCLE_PCA")
scrna <- CellCycleScoring(scrna, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
scrna$G1.Score = 1 - scrna$S.Score - scrna$G2M.Score
scrna$CC.Difference <- scrna$S.Score - scrna$G2M.Score

scrna <- ScaleData(scrna,
                   vars.to.regress = c("nCount_RNA","percent.mt","G2M.Score","S.Score"),
                   features = rownames(scrna))
scrna <- FindVariableFeatures(scrna)
scrna <- RunPCA(scrna, features = VariableFeatures(scrna), nfeatures.print = 10, reduction.name="RegressOut_PCA")

DefaultAssay(scrna) <- "RNA"
scrna <- harmony::RunHarmony(scrna,  "Sample", plot_convergence = "True", reduction = "RegressOut_PCA")
saveRDS(scrna, "save/scrna_dm.rds")
