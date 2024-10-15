suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
mute = function(x) suppressWarnings(suppressMessages(x))

ids =c("353_df",
       "353_si",
       "353_d12",
       "353_d18")


ids2mannual = list(
 "353_df.1" =  "0018183_R1a.csv",
 "353_df.2" =  "0018183_R1b.csv",
 "353_si"   =  "0018183_R3.csv",
 "353_d12"  =  "0018195_R1.csv",
 "353_d18.1" = "0018195_R2a.csv",
 "353_d18.2" = "0018195_R2b.csv"
)




flist = list(
  f353_df  = "data/output-XETG00229__0018183__Region_1__7+18_ips353_dharmafect_ctrl",
  f353_si  = "data/output-XETG00229__0018183__Region_3__7+18_ips353_siRNA_multiplex",
  f353_d12 = "data/output-XETG00229__0018195__Region_1__7+12_ips353",
  f353_d18 = "data/output-XETG00229__0018195__Region_2__7+18_ips353",
)

if(F){## loading each sample to seurat object
  for(x in ids){
    message(date(), " ", x)
    f = flist[[paste0("f", x)]]
    obj <- mute(LoadXenium(f, fov = "fov"))
    saveRDS(obj, glue::glue("save/iPS{x}_obj.Rds"))
  }
}

for(x in ids){## filtering by manually selected region.
  message(date(), " ", x)
  f = flist[[paste0("f", x)]]
  obj <- readRDS(glue::glue("save/iPS{x}_obj.Rds"))
  regs <-names(ids2mannual)[ which(grepl(paste0("^", x), names(ids2mannual)))]
  cell_ids <- c()
  for(reg in regs){
    message(date(), " ", reg)
    df_f = file.path("data/manually_selected_region", ids2mannual[[reg]])
    df = read.csv(df_f, skip=2, stringsAsFactors = F)
    #df$cell_id = paste0(df$cell_id, "_", reg)
    assertthat::assert_that(all(df$Cell.ID%in% rownames(obj@meta.data)))
    cell_ids = c(cell_ids, df$Cell.ID)
  }
  obj <- subset(obj, cells = cell_ids)
  saveRDS(obj, glue::glue("save/iPS{x}_obj_cleaned.Rds"))
}
