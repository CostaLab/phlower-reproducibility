library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)

findids <- function(coordx, coordy, xrange, yrange){
    library(spatstat.utils)
    ids <- which(inside.range(coordx, xrange) & inside.range(coordy, yrange))
    return(ids)
}

ids =c("353_df",
      "353_si",
      "353_d12",
      "353_d18")


flist = list(
  f353_df  = "../data/output-XETG00229__0018183__Region_1__7+18_ips353_dharmafect_ctrl",
  f353_si  = "../data/output-XETG00229__0018183__Region_3__7+18_ips353_siRNA_multiplex",
  f353_d12 = "../data/output-XETG00229__0018195__Region_1__7+12_ips353",
  f353_d18 = "../data/output-XETG00229__0018195__Region_2__7+18_ips353"
)

annotation_353 <- c("0" = "Tubular Prog.",
                    "1" = "Podocytes Prog.",
                    "2" = "Stromal",
                    "3" = "Podocytes",
                    "4" = "Epithelial",
                    "5" = "Neuronal",
                    "6" = "Tubular Prog.",
                    "7" = "Tubular",
                    "8" = "Mesoderm",
                    "9" = "Tubular Prog.",
                    "10" = "Epithelial",
                    "11" = "Stromal Prog.",
                    "12" = "Epithelial",
                    "13" = "Epithelial",
                    "14"= "Tubular Prog.")

color_353 <- c(
  "Podocytes" = "#16FF32",
  "Podocytes Prog." = "#1CBE4F",
  "Stromal" = "#F6222E",
  "Stromal Prog." = "#822E1C",
  "Tubular" = "#3283FE",
  "Tubular Prog." = "#3B00FB",
  "Epithelial" = "#1A1571",
  "Mesoderm" = "#FEAF16",
  "Neuronal" = "#5A5156"
)

#adjust from ggsci::pal_simpsons
sample_color = c(
  iPS353_d12 = "#FD7446FF",
  d18_pos_top = "#FED439FF",
  d18_pos_down = "#675617FF",
  df_pos_top = "#8A9197FF" ,
  df_pos_down = "#3D4043FF",
  iPS353_si = "#D5E4A2FF"
)

annotation_all <- c("0" = "Podocytes",
                    "1" = "Stromal Prog.",
                    "2" = "Tubular Prog.",
                    "3" = "Tubular Prog.",
                    "4" = "Neuronal",
                    "5" = "Podocytes Prog.",
                    "6" = "Tubular",
                    "7" = "Podocytes Prog.",
                    "8" = "Tubular Prog.",
                    "9" = "Mesoderm",
                    "10" = "Stromal",
                    "11" = "Podocytes Prog.",
                    "12" = "Neuronal Prog.",
                    "13" = "Tubular Prog.",
                    "14" = "Stromal", #?
                    "15" = "Stromal", #?
                    "16" = "Mesoderm", #?
                    "17"= "Tubular")



color_all <- c(
  "Podocytes" = "#16FF32",
  "Podocytes Prog." = "#1CBE4F",
  "Stromal" = "#F6222E",
  "Stromal Prog." = "#822E1C",
  "Tubular" = "#3283FE",
  "Tubular Prog." = "#3B00FB",
  "Mesoderm" = "#FEAF16",
  "Neuronal" = "#5A5156",
  "Neuronal Prog." = "#3F161B"#,
)


color_353_highlight_epi <- function(color_353, highlight_celltype="Podocytes"){
  for (celltype in names(color_353)){
    if (stringr::str_detect(celltype, highlight_celltype)){
      #color_353[celltype] <- "red"
      if (highlight_celltype %in% c("Mesoderm", "Neuronal", "Epithelial")){
        color_353[celltype] <- "red"
      } else {
        next
      }
    }else{
      color_353[celltype] <- "lightgrey"
    }
  }
  #if(highlight_celltype == "Epithelial"){
  #  color_353[""] <- "red"
  #}

  return(color_353)
}


color_353_highlight <- function(color_353, highlight_celltype="Podocytes"){
  for (celltype in names(color_353)){
    if (stringr::str_detect(celltype, highlight_celltype)){
      #color_353[celltype] <- "red"
      if (highlight_celltype %in% c("Mesoderm", "Neuronal", "Epithelial")){
        color_353[celltype] <- "red"
      } else {
        next
      }
    }else{
      color_353[celltype] <- "lightgrey"
    }
  }
  return(color_353)
}

scProportionPlotS <- function(project,
                              is_barplot=TRUE,
                              clusterName="Cluster_0.5",
                              condition="Sample",
                              pair=c("a", "b"),
                              color = color_353
){

#TODO:
  #add checking the input parameters's validation

  prop_test <- sc_utils(project)
  #prop_test@meta_data$sample_identity = prop_test@meta_data[, condition]
  #prop_test@meta_data$cluster_identity = prop_test@meta_data[, clusterName]

  prop_test@meta_data <-  prop_test@meta_data[eval(as.name(condition)) %in% pair ]

  prop_test <- permutation_test(prop_test,
                                cluster_identity = clusterName,
                                sample_1 = pair[1],
                                sample_2 = pair[2],
                                sample_identity = condition)

  dff = data.table::as.data.table(prop_test@results$permutation)
  if(is_barplot){
    p <- permutation_barplot(dff, color=color) + ggtitle(glue::glue("{pair[1]} vs {pair[2]}, positive means more {pair[2]}"))
  }else{

    p <- permutation_plot(dff) + ggtitle(glue::glue("{pair[1]} vs {pair[2]}, positive means more {pair[2]}"))
  }
  p
}


PropBarPlotS <- function(project, Cluster="Cluster", condition="Sample"){

  dfm = as.data.frame.table(table(as.vector(project@meta.data[, condition]), as.vector(project@meta.data[, Cluster]))) %>%
                                                dplyr::mutate(proportion=Freq/sum(Freq)) %>%
                                                dplyr::select(-c(Freq))
  colnames(dfm) <- c(condition, Cluster, "proportion")

  dfm <- ddply(dfm, condition, transform,Share=proportion)
  #print(dfm)

  ggplot(dfm) +
          aes(x=!!sym(Cluster), y=proportion,  fill = !!sym(condition), label=round(proportion, 2)) +
          geom_bar(position = "fill", stat = "identity") +
          #geom_text(aes(label = scales::percent(proportion)),position="stack",vjust=+2.1,size=3) +
          scale_y_continuous(labels=scales::percent)
}


PiePlot <- function(project, Cluster, condition=NULL, cols=ggsci::pal_igv()(51), round_n=2, label_size=3, ncol=3){

  ## colors would be wrong, when the cells not present in some clusters
  if (is.null(condition)){
      data = as.data.frame.table(table(project@meta.data[, Cluster])) %>% dplyr::mutate(proportion=round(100.0*Freq/sum(Freq), round_n))
      colnames(data) <- c(Cluster, "cells", "proportion")

      df2 <- data %>%   mutate(csum = rev(cumsum(rev(proportion))),
                                   pos = proportion/2 + lead(csum, 1),
                                   pos = if_else(is.na(pos), proportion/2, pos))

      p <- ggplot(data, aes(x="", y=proportion, fill=!!sym(Cluster))) +
                  geom_bar(stat="identity", width=1, color="white") +
                  coord_polar("y", start=0) + theme_void() +
                  scale_fill_manual(values=cols) +
                  geom_label_repel(data = df2,
                      aes(y = pos, label = paste0(proportion, "%")),
                      size = label_size, nudge_x = 1, show.legend = FALSE)

  }else{
      plist <- list()
      conditions = unique(project@meta.data[, condition])
      for(i in seq_along(conditions)){
        cond = conditions[i]
        meta = project@meta.data
        idx = which(as.vector(meta[, condition]) == cond)
        data = as.data.frame.table(table(meta[idx, Cluster])) %>% dplyr::mutate(proportion=round(100.0*Freq/sum(Freq), round_n))
        colnames(data) <- c(Cluster, "cells", "proportion")
        df2 <- data %>%   mutate(csum = rev(cumsum(rev(proportion))),
                                   pos = proportion/2 + lead(csum, 1),
                                   pos = if_else(is.na(pos), proportion/2, pos))

        px <- ggplot(data, aes(x="", y=proportion, fill=!!sym(Cluster))) +
                  geom_bar(stat="identity", width=1, color="white") +
                  coord_polar("y", start=0) + theme_void() +
                  scale_fill_manual(values=cols) +
                  geom_label_repel(data = df2,
                      aes(y = pos, label = paste0(proportion, "%")),
                      size = label_size, nudge_x = 1, show.legend = FALSE) + ggtitle(cond)

        plist[[cond]] <- px+theme(legend.position='none')
      }

      pl <- ggplot(data, aes(x="", y=proportion, fill=!!sym(Cluster))) +
                  geom_bar(stat="identity", width=1, color="white") +
                  coord_polar("y", start=0) + theme_void() +
                  scale_fill_manual(values=cols) +
                  geom_label_repel(data = df2,
                      aes(y = pos, label = paste0(proportion, "%")),
                      size = label_size, nudge_x = 1, show.legend = FALSE) + theme(legend.direction='vertical')
      leg <- ggpubr::get_legend(pl)

      p <- cowplot::plot_grid(plotlist=plist, ncol=ncol)
      p <- p + ggpubr::as_ggplot(leg) + patchwork::plot_layout(widths = c(7, 1))
  }
  p
}



PiePlotS <- function(project, Cluster, condition=NULL, cols=ggsci::pal_igv()(51), round_n=2, ...){

  ## colors would be wrong, when the cells not present in some clusters
  if (is.null(condition)){
      data = as.data.frame.table(table(project@meta.data[, Cluster])) %>% dplyr::mutate(proportion=round(100.0*Freq/sum(Freq), round_n))
      colnames(data) <- c(Cluster, "cells", "proportion")

      if(is.null(names(cols))){
        #do nothing
        ncols = cols
      }else{
        ncols <- cols[unique(data$Cluster)]
      }
      p <- ggpubr::ggpie(data, x="cells", label = paste0(data$proportion, '%'),
                 fill = Cluster,  ...) + scale_fill_manual(values=ncols) + theme(legend.direction='vertical', legend.position='right')

  }else{
      plist <- list()
      conditions = unique(project@meta.data[, condition])
      for(i in seq_along(conditions)){
        cond = conditions[i]
        meta = project@meta.data
        idx = which(as.vector(meta[, condition]) == cond)
        data = as.data.frame.table(table(meta[idx, Cluster])) %>% dplyr::mutate(proportion=round(100.0*Freq/sum(Freq), round_n))

        if(is.null(names(cols))){
          #do nothing
          ncols = cols
        }else{
          ncols <- cols[unique(data$Cluster)]
        }
        colnames(data) <- c(Cluster, "cells", "proportion")
        px <-ggpubr::ggpie(data, x="cells", label = paste0(data$proportion, '%'),
                   fill = Cluster, ...) + scale_fill_manual(values=ncols)   + ggtitle(cond)

        plist[[cond]] <- px+theme(legend.position='none')
      }

      if(is.null(names(cols))){
        #do nothing
        ncols = cols
      }else{
        ncols <- cols[unique(data$Cluster)]
      }
      pl =ggpubr::ggpie(data, x="cells", fill=Cluster)+ scale_fill_manual(values=ncols)  + theme(legend.direction='vertical')
      leg <- ggpubr::get_legend(pl)

      p <- cowplot::plot_grid(plotlist=plist, ncol=3)
      p <- p + ggpubr::as_ggplot(leg) + patchwork::plot_layout(widths = c(7, 1))
  }
  p
}

