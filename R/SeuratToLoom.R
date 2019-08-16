# process seurat to loom
library(loomR)
library(hdf5r)
library(Seurat)
seur <- readRDS("/Users/jsjoyal/Desktop/SCAP/test_data/seurat_stim_vs_contol/seurat_obj.rds")

assays <- names(seur@assays)

flag <- FALSE
if(any(assays == 'integrated')){
  assays <- assays[-grep('integrated|SCT',assays)]
}else if(any(assays=='SCT')){
    assays <- assays[-grep('integrated|RNA',assays)]
    flag <- TRUE
}

project_dir <- "/Users/jsjoyal/Desktop/SCAP/test_data/seurat_project/"

for(assay in assays){
  DefaultAssay(seur) <- assay
  if(flag == TRUE && assay == "SCT"){
    pseudo.assay <- "RNA"
  }else{
    pseudo.assay <- assay
  }
  filename <- paste0(project_dir,pseudo.assay,".loom")
  loomR::create(filename = filename, data = seur[[assay]]@data, calc.count = F, overwrite = T)
  data <- loomR::connect(filename = filename, mode = "r+")
  data$link_delete('row_attrs/Gene')
  # add metadata'
  meta.data <- seur@meta.data
  colnames(meta.data) <- paste0(colnames(meta.data),"_meta_data")
  if(length(assays)>1){
    remove.assays <- paste(c(paste0("_",assays[-which(assays == pseudo.assay)]),paste0(assays[-which(assays == pseudo.assay)],"_")),collapse = "|")
    rm <- grep(remove.assays,colnames(meta.data), ignore.case = TRUE)
    if(!identical(rm, integer(0))){
      data$add.col.attribute(as.list(meta.data[,-rm]))
    }else{
      data$add.col.attribute(as.list(meta.data))
    }
  }else{
    data$add.col.attribute(as.list(meta.data))
  }
  data$add.row.attribute(list(features = rownames(seur[[pseudo.assay]])))
  # add reductions
  reductions <- tolower(paste0(c("pca","tsne","umap","diff"),"_",pseudo.assay,"_",c(2,2,2,2,3,3,3,3),"d"))
  for(reduc in reductions){
    reduc.data <- seur@reductions[[reduc]]
    if(is.null(reduc.data)) next
    red.df <- as.data.frame(reduc.data@cell.embeddings)
    n <- ncol(red.df)
    if(grepl('pca_',reduc)){
      if(grepl('_2d',reduc)){
        n <- 2
      }else if(grepl('_3d',reduc)){
        n <- 3
      }
    }
    colnames(red.df) <- paste0(reduc,'_',1:n,'_reduction')
    data$add.col.attribute(as.list(red.df[,1:n]))
  }
  data$close_all()
}

