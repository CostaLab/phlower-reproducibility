suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
mute = function(x) suppressWarnings(suppressMessages(x))

ids =c("353_df",
       "353_si",
       "353_d12",
       "353_d18")

obj.list <- list()
for (x in ids){
    obj.list[[x]] = readRDS(glue::glue("save/iPS{x}_obj_cleaned.Rds"))
}

if(F){ ## merge objects to one
  options(future.globals.maxSize = 8000 * 1024^2)
  object <- merge(obj.list[[1]], y = obj.list[ids[-1]], add.cell.ids = paste0("iPS", ids), project = "xenium", merge.data = TRUE)

  object@meta.data$sample = sapply(stringr::str_split(colnames(object), "_"), function(x) paste0(x[1], "_", x[2]))


  object <- saveRDS(object, "save/all_353_cleaned.Rds")

  dev.off()
}

object <- readRDS("save/all_353_cleaned.Rds")
DefaultAssay(object) <- "Xenium"

if(T){
  pdf("viz/all_353_cleaned_VlnPlot_QC_cleaned.pdf", width = 12, height = 8)
  p <- VlnPlot(object, features = c("nFeature_Xenium", "nCount_Xenium"), group.by='sample', ncol = 1, pt.size = 0)
  print(p)
  dev.off()

  options(repr.plot.width=36, repl.plot.height=6)
  genes = c("GRHL2", "WT1", "TWIST1", "PAX3")
  pdf("viz/all_353_cleaned_four_genes.pdf", width = 12, height = 8)
  for (nm in names(obj.list)){
      p <- ImageDimPlot(obj.list[[nm]], fov = "fov", molecules = genes[1:4], nmols = 20000)
      print(p+ggtitle(nm))
  }
  dev.off()
}

Idents(object) <- "Xenium"
object <- mute(subset(object, subset = nCount_Xenium > 0))

#object <- NormalizeData(object)
VariableFeatures(object) <- rownames(object)
object <- SCTransform(object, assay='Xenium')


object <- mute(RunPCA(object))

#object = mute(harmony::RunHarmony(object,  "sample", reduction.use = "pca"))

object <- mute(FindNeighbors(object, reduction = "pca", dims = 1:50))
object <- mute(FindClusters(object, resolution = 0.3))

object <- mute(RunUMAP(object, dims = 1:50, reduction = 'pca'))

options(repr.plot.width=12, repl.plot.height=8)
DimPlot(object, group.by='seurat_clusters', label=T)

options(repr.plot.width=12, repl.plot.height=8)
DimPlot(object, group.by='sample') + ggsci::scale_color_igv()

object$ips <- object$sample
object$ips[object$ips %in% c("12", "15", "18", "df", "nt", "si")] <- "ips353"
#object$ips[object$ips %in% c("nall", "ndf", "nnt", "nsi")] <- "ips15"


DimPlot(object, group.by='ips') + ggsci::scale_color_igv()

#pt, podo, stromal1, stromal234, neuron12, neuron-prog
markers <- c(genes <- c("ESRP1", "BLNK", "DCDC2", "MACC1", "PAX2", "LHX1", "RASSF6", "EMX2", "CACNA2D3", "GRHL2", 
                        "NPHS2", "CLIC5", "ZBTB7C", "WT1", "NPHS1", "PODXL", "NRG3", "PTPRO", "MAFB", "TESC", 
                        "GPC5", "L1TD1", "TWIST1", "SH2D3C", "EDNRA", "LAMA4", "WNT5B", 
                        "SHOX2", "RSPO2", "BEND4", "MKX", "LHX9", "PTGFR", "DLX2", "GABRR1", "FGF14", "CDH12", 
                        "LMX1A", "WSCD2", "RFX4", "WNT10B", "RAB3C", "FYB2", "TMEM255A", "CR1L", "EDAR", "DCC", 
                        "SMPX", "IGSF21", "HTR2C", "PAX7", "MYO18B", "FRMPD1", "MYF6", "MYL1", "CELF2", "MRLN"))

library(ggplot2)
options(repr.plot.width=30, repl.plot.height=5)
DotPlot(object, features = markers, group.by='seurat_clusters') + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))



#rm(obj1, obj2, obj3,obj4,obj5,obj6)
#gc()

saveRDS(object, "save/all_353_pca_cleaned.Rds")



