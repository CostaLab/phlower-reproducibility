library(Seurat)
library(dplyr)
library(DropletUtils)


object <- readRDS(file="save/multiome_mojitoo.rds")

mojitoo <- Embeddings(object[["MOJITOO"]])
write.csv(mojitoo, "save/mojitoo.csv", quote=F)

meta <- object@meta.data
write.csv(meta, "save/meta.csv", quote=F)

rna_ <- GetAssayData(object, assay="RNA", slot="counts")
atac_ <- GetAssayData(object, assay="Peaks", slot="counts")


DropletUtils::write10xCounts(path="save/rna_counts.h5", x=rna_, type="HDF5")
DropletUtils::write10xCounts(path="save/atac_counts.h5", x=atac_, type="HDF5")
