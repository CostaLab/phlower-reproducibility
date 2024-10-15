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

source("util/util.R")
source("util/scProportion.R")

mute = function(x) suppressWarnings(suppressMessages(x))

message(date(), " Loading data")
object <- readRDS("save/all_353_pca_cleaned.Rds")
message(date(), " Data loaded")


#100 genes
genes <- c("ESRP1", "BLNK", "DCDC2", "MACC1", "PAX2", "LHX1", "RASSF6",
           "EMX2", "CACNA2D3", "GRHL2", "NPHS2", "CLIC5", "ZBTB7C", "WT1",
           "NPHS1", "PODXL", "NRG3", "PTPRO", "MAFB", "TESC", "GPC5", "L1TD1",
           "TWIST1", "SH2D3C", "EDNRA", "LAMA4", "WNT5B", "SHOX2", "RSPO2",
            "BEND4", "MKX", "LHX9", "PTGFR", "DLX2", "GABRR1", "FGF14",
            "CDH12", "LMX1A", "WSCD2", "RFX4", "WNT10B", "RAB3C", "FYB2",
            "TMEM255A", "CR1L", "EDAR", "DCC", "SMPX", "IGSF21", "HTR2C",
            "PAX7", "MYO18B", "FRMPD1", "MYF6", "MYL1", "CELF2", "MRLN",
            "FOXC2", "FOXP1", "FOXO3", "MAF", "MAFG", "ERF", "CLOCK", "ETV5",
            "MAX", "FLI1", "NFIA", "NFAT5", "RUNX1", "NFATC4", "NFIB",
            "PAX3", "HOXB8", "LEF1", "HOXC10", "MYOG", "MSC", "MYF5",
            "TBXT", "KDR", "PAPPA2", "SLC12A1", "MSX1", "MAP2", "PDGFRB",
            "COL1A2", "HNF1B", "ZIC2", "RUNX2", "ERBB4", "CPEB1", "CHRM3",
            "WNT1", "POU4F1", "POU3F3", "FRY", "FRMPD4", "RGS7", "CDH6")


if(F){
  library(tidyverse)

  #dff <- read.csv("save/all_meta.csv", row.names=1)
  dff <- object@meta.data

  gdf <- dff %>% group_by(sample) %>% summarise(nCount_Xenium=sum(nCount_Xenium),
                                                nFeature_Xenium=sum(nFeature_Xenium))

  df2 <- gdf %>%
    mutate(csum = rev(cumsum(rev(nCount_Xenium))),
           pos = nCount_Xenium/2 + lead(csum, 1),
           pos = if_else(is.na(pos), nCount_Xenium/2, pos))

  p1 <- ggplot(gdf, aes(x = "", y = nCount_Xenium, fill = fct_inorder(sample))) +
    geom_col(width = 1, color = 1) +
    geom_text(aes(label = nCount_Xenium), position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    guides(fill = guide_legend(title = "sample")) +
    scale_y_continuous(breaks = df2$pos, labels = gdf$sample) +
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 15),
          legend.position = "none", # Removes the legend
          panel.background = element_rect(fill = "white"))



  gdf <- dff %>% group_by(sample) %>% summarise(nFeature_Xenium=sum(nFeature_Xenium),
                                                nFeature_Xenium=sum(nFeature_Xenium))


  df2 <- gdf %>%
    mutate(csum = rev(cumsum(rev(nFeature_Xenium))),
           pos = nFeature_Xenium/2 + lead(csum, 1),
           pos = if_else(is.na(pos), nFeature_Xenium/2, pos))

  p2 <- ggplot(gdf, aes(x = "", y = nFeature_Xenium, fill = fct_inorder(sample))) +
    geom_col(width = 1, color = 1) +
    geom_text(aes(label = nFeature_Xenium), position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    guides(fill = guide_legend(title = "sample")) +
    scale_y_continuous(breaks = df2$pos, labels = gdf$sample) +
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 15),
          legend.position = "none", # Removes the legend
          panel.background = element_rect(fill = "white"))



  gdf <- dff %>% group_by(sample) %>% summarise(n=n())

  df2 <- gdf %>%
    mutate(csum = rev(cumsum(rev(n))),
           pos = n/2 + lead(csum, 1),
           pos = if_else(is.na(pos), n/2, pos))

  p3 <- ggplot(gdf, aes(x = "", y = n, fill = fct_inorder(sample))) +
    geom_col(width = 1, color = 1) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    guides(fill = guide_legend(title = "sample")) +
    scale_y_continuous(breaks = df2$pos, labels = gdf$sample) +
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 15),
          legend.position = "none", # Removes the legend
          panel.background = element_rect(fill = "white"))





  pdf("viz/all_353_cleaned_seurat_clusters/counts.pdf")
  print(p1 + ggtitle("Number of counts"))
  print(p2 + ggtitle("Number of features"))
  print(p3 + ggtitle("Number of cells"))
  dev.off()
}


if(F){
  pdf("viz/all_353_cleaned_seurat_clusters/umap_seurat_clusters.pdf", width=10, height=8)
  p <- DimPlot(object, group.by='sample', label=F) + ggsci::scale_color_igv()
  print(p)
  object$seurat_clusters <- as.character(object$seurat_clusters)
  p <- DimPlot(object, group.by='seurat_clusters', label=T) + ggsci::scale_color_igv()
  print(p)
  dev.off()
}
if(F){
  markers <- c("TBXT","KDR",      #Mesoderm
               "PODXL","NPHS2",   #Podo
               "PAPPA2","SLC12A1",#Tubular
               "MAP2","MSX1",     #Neuron
               "PDGFRB", "COL1A2")#Stromal

  pdf("viz/all_353_cleaned_seurat_clusters/all_dotplot_seurat_clusters.pdf", width=10, height=8)
  p <- DotPlot(object, features = markers, group.by='seurat_clusters') +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  print(p)
  dev.off()
}


if(F){
  pdf("viz/all_353_cleaned_seurat_clusters/ko_genes_barplot.pdf")
  sct_mtx <- GetAssayData(object, assay='SCT', layer='scale.data')
  meta <- object@meta.data

  ko_genes <- c("RFX4", "ZIC2", "PAX3")
  for (ko in ko_genes){
      sub_df<- as.data.frame.matrix(t(sct_mtx[ko,,drop=F]))
      sub_df <- cbind(sub_df, meta[, "sample"])
      mdf <- melt(sub_df)
      colnames(mdf) <- c("sample", "gene", "value")
      mmdf <- mdf %>% group_by(sample) %>%  dplyr::summarize(Mean = mean(value, na.rm=TRUE))
      p <- ggplot(mmdf, aes(x=sample, y=Mean, fill=sample)) +
                          geom_bar(stat="identity") +
                          ggtitle(ko)
      print(p)
  }

  ## exclude d5
  #ko_genes <- c("RFX4", "ZIC2", "PAX3")
  #for (ko in ko_genes){
  #    sub_df<- as.data.frame.matrix(t(sct_mtx[ko,,drop=F]))
  #    sub_df <- cbind(sub_df, meta[, "sample"])
  #    #remove sample d5
  #    #sub_df <- sub_df %>% filter(sample != "d5")
  #    mdf <- melt(sub_df)
  #    colnames(mdf) <- c("sample", "gene", "value")
  #    mmdf <- mdf %>% group_by(sample) %>%  dplyr::summarize(Mean = mean(value, na.rm=TRUE))
  #    p <- ggplot(mmdf, aes(x=sample, y=Mean, fill=sample)) +
  #                        geom_bar(stat="identity") +
  #                        ggtitle(ko)
  #    print(p)
  #}
  #dev.off()



  dev.off()
}

if(F){
  #pt, podo, stromal1, stromal234, neuron12, neuron-prog
  markers_more <- c("ESRP1", "BLNK", "DCDC2", "MACC1", "PAX2", "LHX1", "RASSF6", "EMX2", "CACNA2D3", "GRHL2", #PT
                    "NPHS2", "CLIC5", "ZBTB7C", "WT1", "NPHS1", "PODXL", "NRG3", "PTPRO", "MAFB", "TESC", #PODO
                    "GPC5", "L1TD1", "TWIST1", "SH2D3C", "EDNRA", "LAMA4", "WNT5B",                        #STROMAL1
                    "SHOX2", "RSPO2", "BEND4", "MKX", "LHX9", "PTGFR", "DLX2", "GABRR1", "FGF14", "CDH12", #STROMAL234
                    "LMX1A", "WSCD2", "RFX4", "WNT10B", "RAB3C", "FYB2", "TMEM255A", "CR1L", "EDAR", "DCC", #NEURON12
                    "SMPX", "IGSF21", "HTR2C", "PAX7", "MYO18B", "FRMPD1", "MYF6", "MYL1", "CELF2", "MRLN") #NEURON-PROG

  pdf("viz/all_353_cleaned_seurat_clusters/all_dotplot_seurat_clusters_more.pdf", width=30, height=8)
  p <- DotPlot(object, features = markers_more, group.by='seurat_clusters') +
             theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  print(p)
  dev.off()
}

if (T){ ## scProportionPlotS test
  #all_pairs
  all_pairs <- combn(unique(object$sample), 2, simplify =  F)
  pdf("viz/all_353_cleaned_seurat_clusters/all_annotatioin_proption_test.pdf", width=8, height=3)
  for(x in all_pairs){
      message(date(), ": ", x[1], " vs ", x[2])
      p <- scProportionPlotS(object, "seurat_clusters", 'sample', c(x[1], x[2]))
      print(p)
  }
  dev.off()
}


if (F){ #  boxplot
  pdf("viz/all_353_cleaned_seurat_clusters/box_all_genes_seurat_clusters_all_scaledata.pdf", width = 8, height = 6)
  DefaultAssay(object) <- "SCT"
  sct_mtx <- GetAssayData(object, assay='SCT', layer='scale.data')
  meta <- object@meta.data

  #ko_genes <- c("RFX4", "ZIC2", "PAX3")
  ko_genes <- genes
  for (ko in ko_genes){
      sub_df<- as.data.frame.matrix(t(sct_mtx[ko,,drop=F]))
      sub_df <- cbind(sub_df, meta[, "sample"])
      mdf <- melt(sub_df)
      colnames(mdf) <- c("sample", "gene", "value")
      p <- ggplot(mdf, aes(x=sample, y=value, fill=sample)) +
                          geom_boxplot(outlier.shape=NA) + ggtitle(ko)
      ylim1 = boxplot.stats(mdf$value)$stats[c(1, 5)]
      print(p + coord_cartesian(ylim = ylim1*1.05))
  }
  dev.off()
}



if (F){
  #barplot proportion
  pdf("viz/all_353_cleaned_seurat_clusters/all_353_cleaned_seurat_clusters_proportion.pdf")
  p1 <- PropBarPlotS(object, 'seurat_clusters', 'sample')+ ggsci::scale_fill_rickandmorty()
  p2 <- PropBarPlotS(object,  'sample', 'seurat_clusters')+ ggsci::scale_fill_igv() #+ scale_fill_manual(values=color)
  print(p1)
  print(p2)
  dev.off()

  #piechart proportion
  pdf("viz/all_353_cleaned_seurat_clusters/all_353_cleaned_seurat_clusters_piechart.pdf", width=16, height=8)
  p1 <- PiePlotS(object, 'seurat_clusters', 'sample')
  p2 <- PiePlot(object, 'seurat_clusters', 'sample', cols=color, label_size=3)
  print(p1)
  print(p2)
  dev.off()
}


if (F){ # ridges and vlnplot
  pdf("viz/all_353_cleaned_seurat_clusters/ridges_all_annotatioin_genes_all_scaledata.pdf", width = 8, height = 6)
  for (g in genes){

      p <- mute(RidgePlot(object, g, group.by='sample', layer='scale.data', assay='SCT'))
      mute(print(p))

  }
  dev.off()
  pdf("viz/all_353_cleaned_seurat_clusters/violin_all_genes_seurat_clusters_all_scaledata.pdf", width = 8, height = 6)
  DefaultAssay(object) <- "SCT"
  for (g in genes){
      p <- mute(VlnPlot(object, g, group.by='sample', layer='scale.data', assay='SCT', pt.size=0))
      mute(print(p))

  }
  dev.off()
}
if(F){
  DefaultAssay(object) <- "Xenium"
  pdf("viz/all_353_cleaned_seurat_clusters/Vln_QC_all_353_cleaned_seurat_clusters.pdf")
  p <- VlnPlot(object, group.by='seurat_clusters', features=c("nFeature_Xenium", "nCount_Xenium"), pt.size=0)
  print(p)

  p <- VlnPlot(object, group.by='sample', features=c("nFeature_Xenium", "nCount_Xenium"), pt.size=0)
  print(p)

  dev.off()
}

# clusters in spatial, ignore
if(F){
  obj21 <- readRDS("save/183_df.Rds")
  obj22 <- readRDS("save/183_nt.Rds")
  obj23 <- readRDS("save/183_si.Rds")
#  obj24 <- readRDS("save/195_dall.Rds")
  obj241 <- readRDS("save/195_d5.Rds")
  obj242 <- readRDS("save/195_d12.Rds")
  obj243 <- readRDS("save/195_d18.Rds")

  fobj = function(x) switch(x, "df"=obj21,
                               "nt"=obj22,
                               "si"=obj23,
                               "d5"=obj241,
                               "d12"=obj242,
                               "d18"=obj243,
                                NULL)

  pdf("viz/all_353_cleaned_seurat_clusters/all_353_cleaned_seurat_clusters_spatail.pdf", width=10, height=7)
  for(id in ids){
      #id = "df"
      message(date(), ": ", id)
      xenium.obj = fobj(id)
      nms <- paste0(id, '_', colnames(xenium.obj))
      inter_nms <- intersect(nms, colnames(object))
      #rm xxx_ from the inter_name to do the filtering
      trim_nms <- stringr::str_remove(inter_nms, paste0("^", id, '_'))
      xenium.obj <- subset(xenium.obj, cells=trim_nms)
      seurat_clusters <- object@meta.data[inter_nms, "seurat_clusters"]
      xenium.obj$seurat_clusters <- seurat_clusters
      xenium.obj$seurat_clusters <- factor(xenium.obj$seurat_clusters, levels=names(color))
      new_color <- color

      Idents(xenium.obj) <- "seurat_clusters"
      if (id == "df"){
          cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(3600, 9000), coords = "plot")
          xenium.obj[["top"]] <- cropped.coords

          cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(0, 3600), coords = "plot")
          xenium.obj[["down"]] <- cropped.coords

          p1 <- ImageDimPlot(xenium.obj, fov = "top", axes = TRUE, border.color = "white",
                             border.size = 0.1,  coord.fixed = FALSE, nmols = 10000) +
                             scale_fill_manual(values=new_color)
          p2 <- ImageDimPlot(xenium.obj, fov = "down", axes = TRUE, border.color = "white",
                             border.size = 0.1,  coord.fixed = FALSE, nmols = 10000) +
                             scale_fill_manual(values=new_color)
          print(p1)
          print(p2)
      }else if((id == "nt") | (id == "time")){

          cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 6000), y = c(5000, 12000), coords = "plot")
          xenium.obj[["top"]] <- cropped.coords

          cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 6000), y = c(0, 5000), coords = "plot")
          xenium.obj[["down"]] <- cropped.coords

          Idents(xenium.obj) <- "seurat_clusters"
          p1 <- ImageDimPlot(xenium.obj, fov = "top", axes = TRUE, border.color = "white",
                             border.size = 0.1,  coord.fixed = FALSE, nmols = 10000)+
                             scale_fill_manual(values=new_color)
          p2 <- ImageDimPlot(xenium.obj, fov = "down", axes = TRUE, border.color = "white",
                             border.size = 0.1,  coord.fixed = FALSE, nmols = 10000)+
                             scale_fill_manual(values=new_color)
          print(p1)
          print(p2)
      }else{
          p <- ImageDimPlot(xenium.obj, fov = "fov", axes = TRUE, border.color = "white", cols='polychrome',
                            border.size = 0.1,  coord.fixed = FALSE, nmols = 10000) +
                            scale_fill_manual(values=new_color)
          print(p)
      }
  }
  dev.off()
}


if(F){ ## highlight one cell type
  obj21 <- readRDS("save/183_df.Rds")
  obj22 <- readRDS("save/183_nt.Rds")
  obj23 <- readRDS("save/183_si.Rds")
#  obj24 <- readRDS("save/195_dall.Rds")
  obj241 <- readRDS("save/195_d5.Rds")
  obj242 <- readRDS("save/195_d12.Rds")
  obj243 <- readRDS("save/195_d18.Rds")

  fobj = function(x) switch(x, "df"=obj21,
                               "nt"=obj22,
                               "si"=obj23,
                               "d5"=obj241,
                               "d12"=obj242,
                               "d18"=obj243,
                                NULL)

  pdf("viz/all_353_cleaned_seurat_clusters/all_353_cleaned_seurat_clusters_spatial_highlight.pdf", width=10, height=7)
  highlight_list <- c("Podocytes", "Tubular", "Stromal", "Mesoderm", "Neuronal")
  for (highlight_celltype in highlight_list){
    new_color <- color_highlight(color, highlight_celltype)
    for(id in ids){
          #id = "df"
          message(date(), ": ", id)
          xenium.obj = fobj(id)
          nms <- paste0(id, '_', colnames(xenium.obj))
          inter_nms <- intersect(nms, colnames(object))
          #rm xxx_ from the inter_name to do the filtering
          trim_nms <- stringr::str_remove(inter_nms, paste0("^", id, '_'))
          xenium.obj <- subset(xenium.obj, cells=trim_nms)
          seurat_clusters <- object@meta.data[inter_nms, "seurat_clusters"]
          xenium.obj$seurat_clusters <- seurat_clusters
          xenium.obj$seurat_clusters <- factor(xenium.obj$seurat_clusters, levels=names(color))

          Idents(xenium.obj) <- "seurat_clusters"
          if (id == "df"){
              cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(3600, 9000), coords = "plot")
              xenium.obj[["top"]] <- cropped.coords

              cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(0, 3600), coords = "plot")
              xenium.obj[["down"]] <- cropped.coords

              p1 <- ImageDimPlot(xenium.obj, fov = "top", axes = TRUE, border.color = "white",
                                 border.size = 0.1,  coord.fixed = FALSE, nmols = 10000) +
                                 scale_fill_manual(values=new_color)
              p2 <- ImageDimPlot(xenium.obj, fov = "down", axes = TRUE, border.color = "white",
                                 border.size = 0.1,  coord.fixed = FALSE, nmols = 10000) +
                                 scale_fill_manual(values=new_color)
              print(p1)
              print(p2)
          }else if((id == "nt") | (id == "time")){
              cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 6000), y = c(5000, 12000), coords = "plot")
              xenium.obj[["top"]] <- cropped.coords

              cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 6000), y = c(0, 5000), coords = "plot")
              xenium.obj[["down"]] <- cropped.coords

              Idents(xenium.obj) <- "seurat_clusters"
              p1 <- ImageDimPlot(xenium.obj, fov = "top", axes = TRUE, border.color = "white",
                                 border.size = 0.1,  coord.fixed = FALSE, nmols = 10000)+
                                 scale_fill_manual(values=new_color)
              p2 <- ImageDimPlot(xenium.obj, fov = "down", axes = TRUE, border.color = "white",
                                 border.size = 0.1,  coord.fixed = FALSE, nmols = 10000)+
                                 scale_fill_manual(values=new_color)
              print(p1)
              print(p2)
          }else{
              p <- ImageDimPlot(xenium.obj, fov = "fov", axes = TRUE, border.color = "white",
                                border.size = 0.1,  coord.fixed = FALSE, nmols = 10000) +
                                scale_fill_manual(values=new_color)
              print(p)
          }
      }

      }
  dev.off()
}


