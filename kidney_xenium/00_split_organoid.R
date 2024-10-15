library(Seurat)


ids =c("353_df",
       "353_si",
       "353_d12",
       "353_d18")


#cells_of_the_cut = cropped.coords@boundaries$centroids@cells
for (id in ids){
  message(date(), " ", id)

  if (id == "353_df"){
      xenium.obj <- readRDS(glue::glue("save/iPS{id}_obj_cleaned.Rds"))

      cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(5000, 9000), coords = "plot")
      xenium.obj[["top"]] <- cropped.coords
      cells_of_the_cut_1 = paste0("iPS353_df_", cropped.coords@boundaries$centroids@cells)

      cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(0, 5000), coords = "plot")
      xenium.obj[["down"]] <- cropped.coords
      cells_of_the_cut_2 = paste0("iPS353_df_", cropped.coords@boundaries$centroids@cells)

      dff = data.frame(bc=c(cells_of_the_cut_1, cells_of_the_cut_2), pos = c(rep("top", length(cells_of_the_cut_1)), rep("down", length(cells_of_the_cut_2))))
      write.csv(dff, glue::glue("save/iPS{id}_obj_cleaned_split.csv"), row.names = F, quote = F)

  }else if(id == "353_d18"){
      xenium.obj <- readRDS(glue::glue("save/iPS{id}_obj_cleaned.Rds"))

      cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(4000, 9000), coords = "plot")
      xenium.obj[["top"]] <- cropped.coords
      cells_of_the_cut_1 = paste0("iPS353_d18_", cropped.coords@boundaries$centroids@cells)

      cropped.coords <- Crop(xenium.obj[["fov"]], x = c(0, 3600), y = c(0, 4000), coords = "plot")
      xenium.obj[["down"]] <- cropped.coords
      cells_of_the_cut_2 = paste0("iPS353_d18_", cropped.coords@boundaries$centroids@cells)

      dff = data.frame(bc=c(cells_of_the_cut_1, cells_of_the_cut_2), pos = c(rep("top", length(cells_of_the_cut_1)), rep("down", length(cells_of_the_cut_2))))
      write.csv(dff, glue::glue("save/iPS{id}_obj_cleaned_split.csv"), row.names = F, quote = F)


  }else{
      message(id, " doesn't need split")
  }
}






