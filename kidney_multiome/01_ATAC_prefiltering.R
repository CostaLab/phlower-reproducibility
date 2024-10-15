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

## set parameters
set.seed(42)
addArchRThreads(threads = 100)
addArchRGenome("hg38")


## Creating Arrow Files

minTSS <- 6
minFrags <- 2500

for(i in 1:4){
   inputFiles <- glue("data/S{i}/outs/atac_fragments.tsv.gz")
   names(inputFiles) <- glue("S{i}")

   ArrowFiles <- createArrowFiles(
            inputFiles = inputFiles,
            sampleNames = names(inputFiles),
            outputNames = names(inputFiles),
            minTSS = minTSS,
            minFrags = minFrags,
            maxFrags = 1e+05,
            QCDir = "QualityControl",
            addTileMat = TRUE,
            addGeneScoreMat = TRUE)
}

ArrowFiles <- c("S1.arrow",
                "S2.arrow",
                "S3.arrow",
                "S4.arrow")


proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "organoid",
  showLogo = FALSE,
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)


dir.create("save")
meta <- proj@cellColData
write.table(meta, "save/cellColData.tsv", sep="\t")
saveArchRProject(ArchRProj = proj, 'save/ATAC',load = FALSE)

