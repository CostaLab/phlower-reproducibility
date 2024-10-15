library(Seurat)
library(WriteXLS)

# Load the Seurat object
object <- readRDS("save/all_353_cleaned_pca_annotation.Rds")

# Find markers
Idents(object) <- "seurat_clusters"

## fix the seurat bug that DefaultAssay is not "RNA"
len = length(object@assays$SCT@SCTModel.list)
for(i in 1:len){
  slot(object = object@assays$SCT@SCTModel.list[[i]], name="umi.assay")<-"Xenium"
}
assertthat::assert_that( all(unlist(SCTResults(object=object, slot="umi.assay")) == "Xenium"))


object <- PrepSCTFindMarkers(object, verbose = FALSE)
markers <- FindAllMarkers(object, only.pos = TRUE, logfc.threshold = 0)

# Save the markers
saveRDS(markers, "save/all_353_cleaned_pca_annotation_markers.Rds")

df.list <- split(markers, markers$cluster)
WriteXLS::WriteXLS(df.list, "save/all_353_cleaned_pca_annotation_markers.xlsx")
