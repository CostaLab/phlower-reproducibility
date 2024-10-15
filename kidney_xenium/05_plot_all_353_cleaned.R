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


object <- readRDS("save/all_353_cleaned_pca_annotation.Rds")


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



if(F){ #  plot showing a celltype each samples

  pdf("viz/all_353_cleaned/sample_violin_in_each_annotation.pdf", width=10, height=10)
  Idents(object) <-"annotation"
  tbl <- table(object$annotation, object$sample)
  dff <-  as.data.frame.table(tbl)
  names(dff) <- c("celltype", "sample", "count")
  for (ct in unique(object$annotation)){
    #p <- VlnPlot(object, group.by= "sample", pt.size=0, idents=ct) + ggtitles(ct)
    p <- ggplot(dff[dff$celltype == ct,], aes(x=sample, y=count)) + geom_violin() + geom_boxplot(width=0.1) + ggtitle(ct)
    print(p)
  }
  dev.off()

}



if(F){
  #library(tidyverse)
  dff <- object@meta.data[, c("nCount_Xenium", "nFeature_Xenium", "sample")]

  gdf <- dff %>% dplyr::group_by(sample) %>% dplyr::summarise(nCount_Xenium=sum(nCount_Xenium),
                                                nFeature_Xenium=sum(nFeature_Xenium))

  df2 <- gdf %>%
    mutate(csum = rev(cumsum(rev(nCount_Xenium))),
           pos = nCount_Xenium/2 + lead(csum, 1),
           pos = if_else(is.na(pos), nCount_Xenium/2, pos))
  message("p1")
  p1 <- ggplot(gdf, aes(x = "", y = nCount_Xenium, fill = forcats::fct_inorder(sample))) +
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



  gdf <- dff %>% dplyr::group_by(sample) %>% dplyr::summarise(nFeature_Xenium=sum(nFeature_Xenium),
                                                nFeature_Xenium=sum(nFeature_Xenium))


  df2 <- gdf %>%
    mutate(csum = rev(cumsum(rev(nFeature_Xenium))),
           pos = nFeature_Xenium/2 + lead(csum, 1),
           pos = if_else(is.na(pos), nFeature_Xenium/2, pos))

  message("p2")
  p2 <- ggplot(gdf, aes(x = "", y = nFeature_Xenium, fill = forcats::fct_inorder(sample))) +
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



  gdf <- dff %>% dplyr::group_by(sample) %>% dplyr::summarise(n=n())

  df2 <- gdf %>%
    mutate(csum = rev(cumsum(rev(n))),
           pos = n/2 + lead(csum, 1),
           pos = if_else(is.na(pos), n/2, pos))

  message("p3")
  p3 <- ggplot(gdf, aes(x = "", y = n, fill = forcats::fct_inorder(sample))) +
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





  pdf("viz/all_353_cleaned/counts.pdf")
  print(p1 + ggtitle("Number of counts"))
  print(p2 + ggtitle("Number of features"))
  print(p3 + ggtitle("Number of cells"))
  dev.off()
}


if(F){
  message(date(), " Plotting umap annotation")
  pdf("viz/all_353_cleaned/umap_annotation.pdf", width=10, height=8)
  p <- DimPlot(object, group.by='annotation', label=T) + scale_color_manual(values = color_353)
  print(p)

  p <- DimPlot(object, group.by='sample', label=F) + ggsci::scale_color_igv()
  print(p)
  dev.off()


  pdf("viz/all_353_cleaned/umap_split_condition", width=50, height=8)
  p <- DimPlot(object, group.by='annotation', split.by='sample', label=T) + ggsci::scale_color_igv()
  print(p)
  dev.off()

  p.list <- list()
  for (i in 1:length(unique(object$sample))){
    p <- DimPlot(object, group.by='annotation',
                 label=T,
                 cells = which(object$sample == unique(object$sample)[i])) +
                        ggsci::scale_color_igv()
    p.list[[i]] <- p + ggtitle(unique(object$sample)[i])
  }

  pdf("viz/all_353_cleaned/umap_annotation_split_sample.pdf", width=13, height=20)
  print(cowplot::plot_grid(plotlist = p.list, ncol = 2))
  dev.off()
}


if(F){
  markers <- c("TBXT","KDR",      #Mesoderm
             "PODXL","NPHS2",   #Podo
             "PAPPA2","SLC12A1",#Tubular
             "MAP2","MSX1",     #Neuron
             "PDGFRB", "COL1A2")#Stromal

  pdf("viz/all_353_cleaned/all_353_cleaned_dotplot_annotation.pdf", width=10, height=8)
  p <- DotPlot(object, features = markers, group.by='annotation') +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  print(p)
  dev.off()
}


if(F){
  message(date(), " Plotting ko_genes_barplot")
  pdf("viz/all_353_cleaned/ko_genes_barplot.pdf")
  sct_mtx <- GetAssayData(object, assay='SCT', layer='scale.data')
  meta <- object@meta.data

  ko_genes <- c("RFX4", "ZIC2", "PAX3")
  for (ko in ko_genes){
      sub_df <- as.data.frame.matrix(t(sct_mtx[ko,,drop=F]))
      sub_df <- cbind(sub_df, meta[, "sample"])
      mdf <- melt(sub_df)
      colnames(mdf) <- c("sample", "gene", "value")
      mmdf <- mdf %>% group_by(sample) %>%  dplyr::summarize(Mean = mean(value, na.rm=TRUE))
      p <- ggplot(mmdf, aes(x=sample, y=Mean, fill=sample)) +
                          geom_bar(stat="identity") +
                          ggtitle(ko)
      print(p+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  }


  dev.off()
}

if(F){
  message(date(), " Plotting markers_more")
  markers_more <- c("ESRP1", "BLNK", "DCDC2", "MACC1", "PAX2", "LHX1", "RASSF6", "EMX2", "CACNA2D3", "GRHL2",
                          "NPHS2", "CLIC5", "ZBTB7C", "WT1", "NPHS1", "PODXL", "NRG3", "PTPRO", "MAFB", "TESC",
                          "GPC5", "L1TD1", "TWIST1", "SH2D3C", "EDNRA", "LAMA4", "WNT5B",
                          "SHOX2", "RSPO2", "BEND4", "MKX", "LHX9", "PTGFR", "DLX2", "GABRR1", "FGF14", "CDH12",
                          "LMX1A", "WSCD2", "RFX4", "WNT10B", "RAB3C", "FYB2", "TMEM255A", "CR1L", "EDAR", "DCC",
                          "SMPX", "IGSF21", "HTR2C", "PAX7", "MYO18B", "FRMPD1", "MYF6", "MYL1", "CELF2", "MRLN")
  pdf("viz/all_353_cleaned/all_353_dotplot_annotation_more.pdf", width=30, height=5)
  p <- DotPlot(object, features = markers_more, group.by='annotation') + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  print(p)
  dev.off()

}

if (T){ ## scProportionPlotS test
  #all_pairs
  message(date(), " Plotting scProportionPlotS")
  all_pairs <- combn(unique(object$sample), 2, simplify =  F)
  all_pairs <- list(
    c("iPS353_df", "iPS353_si"),
    c("iPS353_d18", "iPS353_si")
  )
  pdf("viz/all_353_cleaned/all_353_cleaned_annotation_proption_test.pdf", width=8, height=3)
  for(x in all_pairs){
      message(date(), ": ", x[1], " vs ", x[2])
      p <- scProportionPlotS(object, "annotation", 'sample', c(x[1], x[2]))
      print(p)
  }
  dev.off()
}


if (F){ #  boxplot
  message(date(), " Plotting boxplot")
  pdf("viz/all_353_cleaned/box_all_353_cleaned_genes_annotation_all_scaledata.pdf", width = 8, height = 6)
  DefaultAssay(object) <- "SCT"
  sct_mtx <- GetAssayData(object, assay='SCT', layer='scale.data')
  meta <- object@meta.data

  ko_genes <- c("RFX4", "ZIC2", "PAX3")
  #ko_genes <- genes
  for (ko in ko_genes){
      sub_df <- as.data.frame.matrix(t(sct_mtx[ko,,drop=F]))
      sub_df <- cbind(sub_df, meta[, "sample"])
      names(sub_df) <- c(ko, 'sample')
      mdf <- melt(sub_df)
      colnames(mdf) <- c("sample", "gene", "value")
      p <- ggplot(mdf, aes(x=sample, y=value, fill=sample)) +
                          geom_boxplot(outlier.shape=NA) +                                                                                                         ggtitle(ko)
      ylim1 = boxplot.stats(mdf$value)$stats[c(1, 5)]
      p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      print(p + coord_cartesian(ylim = ylim1*1.05))
  }
  dev.off()
}

if (F){
  #barplot proportion
  message(date(), " Plotting barplot proportion")
  pdf("viz/all_353_cleaned/all_353_cleaned_annotation_proportion.pdf")
  p1 <- PropBarPlotS(object, 'annotation', 'sample')+ ggsci::scale_fill_rickandmorty()
  p2 <- PropBarPlotS(object,  'sample', 'annotation')+ ggsci::scale_fill_igv() + scale_fill_manual(values=color_353)
  print(p1+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  print(p2+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()

  #piechart proportion
  pdf("viz/all_353_cleaned/all_353_cleaned_annotation_piechart.pdf", width=16, height=8)
  p1 <- PiePlotS(object, 'annotation', 'sample')
  p2 <- PiePlot(object, 'annotation', 'sample', cols=color_353, label_size=3)
  print(p1)
  print(p2)
  dev.off()
}


if (F){ # ridges and vlnplot
  message(date(), " Plotting ridges and vlnplot")
  pdf("viz/all_353_cleaned/ridges_all_353_cleaned_annotation_genes_all_scaledata.pdf", width = 8, height = 6)
  for (g in genes){

      p <- mute(RidgePlot(object, g, group.by='sample', layer='scale.data', assay='SCT'))
      mute(print(p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))))

  }
  dev.off()
  pdf("viz/all_353_cleaned/violin_all_genes_annotation_all_scaledata.pdf", width = 8, height = 6)
  DefaultAssay(object) <- "SCT"
  for (g in genes){
      p <- mute(VlnPlot(object, g, group.by='sample', layer='scale.data', assay='SCT', pt.size=0))
      mute(print(p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))))

  }
  dev.off()
}
if(F){
  message(date(), " Plotting Vln QC")
  DefaultAssay(object) <- "Xenium"
  pdf("viz/all_353_cleaned/Vln_QC_all_353_cleaned_annotation.pdf")
  p <- VlnPlot(object, group.by='annotation', features=c("nFeature_Xenium", "nCount_Xenium"), pt.size=0)
  print(p +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

  p <- VlnPlot(object, group.by='sample', features=c("nFeature_Xenium", "nCount_Xenium"), pt.size=0)
  print(p +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

  dev.off()
}

if(F){
  message(date(), " Plotting all spatial")
  obj.list <- lapply(ids, function(x) readRDS(glue::glue("save/iPS{x}_obj_cleaned.Rds")))
  names(obj.list) <- ids

  pdf("viz/all_353_cleaned/all_353_cleaned_annotation_spatial.pdf", width=10, height=7)
  for(id in ids){
      #id = "df"
      message(date(), ": ", id)
      xenium.obj = obj.list[[id]]
      nms <- paste0("iPS", id, '_', colnames(xenium.obj))
      inter_nms <- intersect(nms, colnames(object))
      #rm xxx_ from the inter_name to do the filtering
      trim_nms <- stringr::str_remove(inter_nms, paste0("^iPS", id, '_'))
      xenium.obj <- subset(xenium.obj, cells=trim_nms)
      annotation <- object@meta.data[inter_nms, "annotation"]
      xenium.obj$annotation <- annotation
      xenium.obj$annotation <- factor(xenium.obj$annotation, levels=names(color_353))
      new_color <- color_353

      Idents(xenium.obj) <- "annotation"



      if (id == "353_df"){
          cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(5000, 9000), coords = "plot")
          xenium.obj[["top"]] <- cropped.coords

          cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(0, 5000), coords = "plot")
          xenium.obj[["down"]] <- cropped.coords

          p1 <- ImageDimPlot(xenium.obj, fov = "top", axes = TRUE, border.color = "white",
                             border.size = 0.1,  coord.fixed = FALSE, nmols = 10000) +
                             scale_fill_manual(values=new_color)
          p2 <- ImageDimPlot(xenium.obj, fov = "down", axes = TRUE, border.color = "white",
                             border.size = 0.1,  coord.fixed = FALSE, nmols = 10000) +
                             scale_fill_manual(values=new_color)
          print(p1 + ggtitle("iPS353_df 1"))
          print(p2 + ggtitle("iPS353_df 2"))
      }else if(id == "353_d18"){

          cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(4000, 9000), coords = "plot")
          xenium.obj[["top"]] <- cropped.coords

          cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(0, 4000), coords = "plot")
          xenium.obj[["down"]] <- cropped.coords

          Idents(xenium.obj) <- "annotation"
          p1 <- ImageDimPlot(xenium.obj, fov = "top", axes = TRUE, border.color = "white",
                             border.size = 0.1,  coord.fixed = FALSE, nmols = 10000)+
                             scale_fill_manual(values=new_color)
          p2 <- ImageDimPlot(xenium.obj, fov = "down", axes = TRUE, border.color = "white",
                             border.size = 0.1,  coord.fixed = FALSE, nmols = 10000) +
                             scale_fill_manual(values=new_color)
          print(p1+ ggtitle("iPS353_d18 1"))
          print(p2+ ggtitle("iPS353_d18 2"))
        }else{
          p <- ImageDimPlot(xenium.obj, fov = "fov", axes = TRUE, border.color = "white", cols='polychrome',
                            border.size = 0.1,  coord.fixed = FALSE, nmols = 10000) +
                            scale_fill_manual(values=new_color)
          print(p+ggtitle(paste0("iPS", id)))
      }
  }
  dev.off()
}


if(F){## highlight one cell type
  message("highlight one cell type")
  obj.list <- lapply(ids, function(x) readRDS(glue::glue("save/iPS{x}_obj_cleaned.Rds")))
  names(obj.list) <- ids

  pdf("viz/all_353_cleaned/all_353_cleaned_annotation_spatial_highlight.pdf", width=10, height=7)
  highlight_list <- c("Podocytes", "Tubular", "Stromal", "Mesoderm", "Neuronal", "Epithelial")
  for (highlight_celltype in highlight_list){
    new_color <- color_353_highlight(color_353, highlight_celltype)
    for(id in ids){

          message(date(), ": ", id)
          xenium.obj = obj.list[[id]]
          nms <- paste0("iPS", id, '_', colnames(xenium.obj))
          inter_nms <- intersect(nms, colnames(object))
          #rm xxx_ from the inter_name to do the filtering
          trim_nms <- stringr::str_remove(inter_nms, paste0("^iPS", id, '_'))
          xenium.obj <- subset(xenium.obj, cells=trim_nms)
          annotation <- object@meta.data[inter_nms, "annotation"]
          xenium.obj$annotation <- annotation
          xenium.obj$annotation <- factor(xenium.obj$annotation, levels=names(color_353))

          Idents(xenium.obj) <- "annotation"
          if (id == "353_df"){
              cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(5000, 9000), coords = "plot")
              xenium.obj[["top"]] <- cropped.coords

              cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(0, 5000), coords = "plot")
              xenium.obj[["down"]] <- cropped.coords

              p1 <- ImageDimPlot(xenium.obj, fov = "top", axes = TRUE, border.color = "white",
                                 border.size = 0.1,  coord.fixed = FALSE, nmols = 10000) +
                                 scale_fill_manual(values=new_color)
              p2 <- ImageDimPlot(xenium.obj, fov = "down", axes = TRUE, border.color = "white",
                                 border.size = 0.1,  coord.fixed = FALSE, nmols = 10000) +
                                 scale_fill_manual(values=new_color)
              print(p1+ ggtitle("iPS353_df 1"))
              print(p2+ ggtitle("iPS353_df 2"))
          }else if(id == "353_d18"){
              cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(4000, 9000), coords = "plot")
              xenium.obj[["top"]] <- cropped.coords

              cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(0, 4000), coords = "plot")
              xenium.obj[["down"]] <- cropped.coords

              Idents(xenium.obj) <- "annotation"
              p1 <- ImageDimPlot(xenium.obj, fov = "top", axes = TRUE, border.color = "white",
                                 border.size = 0.1,  coord.fixed = FALSE, nmols = 10000)+
                                 scale_fill_manual(values=new_color)
              p2 <- ImageDimPlot(xenium.obj, fov = "down", axes = TRUE, border.color = "white",
                                 border.size = 0.1,  coord.fixed = FALSE, nmols = 10000) +
                                 scale_fill_manual(values=new_color)
              print(p1+ ggtitle("iPS353_d18 1"))
              print(p2+ ggtitle("iPS353_d18 2"))
          }else{
                p <- ImageDimPlot(xenium.obj, fov = "fov", axes = TRUE, border.color = "white", cols='polychrome',
                              border.size = 0.1,  coord.fixed = FALSE, nmols = 10000) +
                              scale_fill_manual(values=new_color)
                print(p+ ggtitle(paste0("iPS", id)))
            }
        }

      }
  dev.off()
}


