# process seurat to loom
library(loomR)
library(hdf5r)
seur <- readRDS("/Users/jsjoyal/Desktop/SCAP/test_data/pbmc.rds")

assays <- names(seur@assays)

project_dir <- "/Users/jsjoyal/Desktop/SCAP/test_data/test_project_2/"

for(assay in assays){
  DefaultAssay(seur) <- assay
  filename <- paste0(project_dir,assay,".loom")
  loomR::create(filename = filename, data = seur[[assay]]@data, calc.count = F, overwrite = T)
  data <- loomR::connect(filename = filename, mode = "r+")
  data$link_delete('row_attrs/Gene')
  # add metadata
  if(length(assay)>1){
    remove.assays <- paste(c(paste0("_",assays[-which(assays == assay)]),paste0(assays[-which(assays == assay)],"_")),collapse = "|")
    data$add.col.attribute(as.list(seur@meta.data[,-c(grep(remove.assays,colnames(seur@meta.data)))]))
  }else{
    data$add.col.attribute(as.list(seur@meta.data))
  }
  data$add.row.attribute(list(features = rownames(seur[[assay]])))
  # add reductions
  reductions <- tolower(paste0(c("pca","tsne","umap","diff"),"_",assay,"_",c(2,2,2,2,3,3,3,3),"d"))
  for(reduc in reductions){
    reduc.data <- seur@reductions[[reduc]]
    if(is.null(reduc.data)) next
    if(grepl('3d',reduc)){
      data$add.col.attribute(as.list(as.data.frame(reduc.data@cell.embeddings)[,1:3]))
    }else{
      data$add.col.attribute(as.list(as.data.frame(reduc.data@cell.embeddings)[,1:2]))
    }
  }
  data$close_all()
}
