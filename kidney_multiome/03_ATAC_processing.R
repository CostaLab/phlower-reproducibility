library(tidyr)
library(dplyr)
library(plotly)
library(cluster)
library(cowplot)
library(gridExtra)
library(viridis)
library(GenomicRanges)
library(GenomeInfoDb)
library(data.table)
library(ArchR)
library(Seurat)
library(glue)
library(BSgenome.Hsapiens.UCSC.hg38)

## set parameters
set.seed(42)
addArchRThreads(threads = 100)
addArchRGenome("hg38")


proj <- loadArchRProject('save/ATAC')
meta <- read.csv("save/cellColData_confiltered.tsv", sep = "\t", row.names = 1)
proj <- proj[rownames(meta), ]
proj <- addIterativeLSI(ArchRProj = proj,
                        name = "IterativeLSI",
                        force = TRUE)

proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)

## cluster for peakset
proj <- addClusters(
    input = proj,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters")
proj <- addReproduciblePeakSet(ArchRProj = proj, groupBy = "Clusters")
proj <- addPeakMatrix(proj)

saveArchRProject(ArchRProj = proj, 'save/ATAC_dm',load = FALSE)
