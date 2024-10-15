mute <- suppressPackageStartupMessages
mute(library(Seurat))
mute(library(ArchR))
mute(library(MOJITOO))
mute(library(tidyverse))
mute(library(BSgenome.Hsapiens.UCSC.hg38))
mute(library(BiocParallel))
mute(library(chromVAR))
mute(library(motifmatchr))
mute(library(stringr))
mute(library(GenomicRanges))
mute(library(IRanges))
mute(library(SummarizedExperiment))
mute(library(JASPAR2020))
mute(library(TFBSTools))
mute(library(Signac))
mute(library(patchwork))
mute(library(dplyr))


set.seed(42)
addArchRThreads(threads = 100)
addArchRGenome("hg38")


addChromVar <- function(object, assay="Peaks",
                                newassay="peaks",
                                sep=c("-", "-"),
                                genome=BSgenome.Hsapiens.UCSC.hg38,
                                genomeName="hg38"){
  # create chromatin assay
  object[[newassay]] <- CreateChromatinAssay(counts=GetAssayData(object[[assay]],slot="counts"),
                                             sep = sep,
                                             genome=genomeName)
  DefaultAssay(object) <- newassay

  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
  )

  object <- AddMotifs(
    object = object,
    genome = genome,
    pfm = pfm
  )


  object <- RunChromVAR(
    object = object,
    genome = genome
  )

  DefaultAssay(object) <- 'chromvar'

  return(object)
}


proj <- loadArchRProject("save/ATAC_dm")
object <- readRDS("save/scrna_dm.rds")
cnms <- colnames(object)
proj <- proj[cnms, ]


ArchR_Harmony <- getDimRed(proj, reduction="Harmony")[cnms, ]
object <- setDimRed(object, mtx=ArchR_Harmony, reduction.name="ArchR_Harmony", key="aharmony")



object <- mojitoo(
     object=object,
     reduction.list = list("harmony", "ArchR_Harmony"),
     dims.list = list(1:50, 1:30),
     is.reduction.center=T,
     reduction.name='MOJITOO',
     assay="RNA"
)

se <- getMatrixFromProject(proj, 'PeakMatrix')
pmtx <- assay(se)
rnms <- GRangesToString(rowRanges(se), sep=c("-", "-"))
pmtx <-pmtx[, colnames(object)]
rownames(pmtx) <- rnms
object[["Peaks"]] <- CreateAssayObject(pmtx)


## add chromvar motif
object <- addChromVar(object)
assertthat::assert_that(all(names(unlist(object@assays$peaks@motifs@motif.names)) == rownames(object)))
object@assays$chromvar@meta.features$tf <- unlist(object@assays$peaks@motifs@motif.names)


## clustering
DefaultAssay(object) <- "RNA"
object <- FindNeighbors(object, reduction="MOJITOO", dims=1:30)
object <- FindClusters(object, resolution=2.5, graph.name="RNA_nn")
object$MOJITOO_cluster <- object$seurat_clusters

saveRDS(object, "save/multiome_mojitoo.rds")
