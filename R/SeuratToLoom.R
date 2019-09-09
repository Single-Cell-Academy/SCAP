# # process seurat to loom
# 
# 
# seuratToLoom <- function(obj, dir){
#   #library(loomR)
#   #library(hdf5r)
#   #library(Seurat)
#   
#   seur <- readRDS(obj)
#   
#   if(!grepl('^seurat',class(seur)[1],ignore.case = T)){
#     showNotification('Error: The selected object is of class ', class(seur)[1], ' but must be of class Seurat', type = 'error')
#   }
#   
#   current_version <- 3
#   if(strsplit(as.character(seur@version), split = '\\.')[[1]][1]<current_version){
#     showNotification(paste0("Warning: The selected seurat object is out of date (version: ", seur@version, "). Updating to object now..."), type = 'error')
#     seur <- UpdateSeuratObject(seur)
#     showNotification(paste0("Update complete"), type = 'message')
#   }
#     
#   project_dir <- paste0(dir,'/')
#   
#   #down_sample <- round(6000/length(unique(seur$seurat_clusters)))
#   #seur <- subset(seur, downsample = down_sample, idents = 'seurat_clusters')
#   # if(dim(seur)[2]>6000){
#   #   seur <- subset(seur, cells = sample(Cells(seur), 6000))
#   # }
#   
#   assays <- names(seur@assays)
#   
#   flag <- FALSE
#   if(any(assays == 'integrated')){
#     assays <- assays[-grep('integrated|SCT',assays)]
#   }else if(any(assays=='SCT')){
#     assays <- assays[-grep('integrated|RNA',assays)]
#     flag <- TRUE
#   }
#   
#   for(assay in assays){
#     DefaultAssay(seur) <- assay
#     if(flag == TRUE && assay == "SCT"){
#       pseudo.assay <- "RNA"
#     }else{
#       pseudo.assay <- assay
#     }
#     filename <- paste0(project_dir,pseudo.assay,".loom")
#     loomR::create(filename = filename, data = seur[[assay]]@data, calc.count = F, overwrite = T)
#     data <- loomR::connect(filename = filename, mode = "r+")
#     data$link_delete('row_attrs/Gene')
#     
#     # add metadata
#     meta.data <- seur@meta.data
#     colnames(meta.data) <- paste0(colnames(meta.data),"_meta_data")
#     
#     data$add.col.attribute(as.list(meta.data))
#     data$add.row.attribute(list(features = rownames(seur[[pseudo.assay]])))
#     
#     # add reduction embeddings
#     reduction_names <- names(seur@reductions)
#     for(i in 1:length(reduction_names)){
#       reductions <- as.data.frame(seur@reductions[[i]]@cell.embeddings)
#       assay_used <- tolower(seur@reductions[[i]]@assay.used)
#       if(!grepl(paste0('(?![a-z])(?<![a-z])(',assay_used,')'),reduction_names[i],perl=T,ignore.case = T)){
#         reduction_names[i] <- paste0(reduction_names[i],'_',assay_used)
#       }
#       n <- ncol(reductions)
#       if(n>3){
#         for(j in 2:3){
#           tmp <- reductions[,1:j]
#           tmp_name <- reduction_names[i]
#           if(!grepl(paste0('(?![a-z])(?<![a-z])(',j,'d)'),tmp_name,perl=T,ignore.case = T)){
#             tmp_name <- paste0(tmp_name,'_',j,'d')
#           }
#           if(grepl(' ', tmp_name)){
#             tmp_name <- gsub(' ','_',tmp_name)
#           }
#           colnames(tmp) <- paste0(tmp_name,'_',1:j,'_reduction')
#           data$add.col.attribute(as.list(tmp))
#         }
#       }else{
#         if(!grepl(paste0('(?![a-z])(?<![a-z])(',n,'d)'),reduction_names[i],perl=T,ignore.case = T)){
#           reduction_names[i] <- paste0(reduction_names[i],'_',n,'d')
#         }
#         if(grepl(' ', reduction_names[i])){
#           reduction_names[i] <- gsub(' ','_',reduction_names[i])
#         }
#         colnames(reductions) <- paste0(reduction_names[i],'_',1:n,'_reduction')
#         data$add.col.attribute(as.list(reductions))
#       }
#     }
#     data$close_all()
#   }
# }
# 
# 
# seurat_reduc_key_name <- function(object, reduction_name = NULL, key_as_name = TRUE, key_name = NULL){
#   if(reduction_name==FALSE) stop('ERROR: you must specify the reduction name')
#   if(key_as_name == FALSE && is.null(key_name)){
#     stop('ERROR: If you do not want to set the key name to the reduction name, a key name must be specified')
#   }else if(key_as_name == TRUE && is.null(key_name)){
#     print('setting key_name to reduction_name...')
#     key_name <- reduction_name
#   }
#   embeddings <- object@reductions[[reduction_name]]@cell.embeddings
#   n <- ncol(embeddings)
#   colnames(embeddings) <- paste0(key_name,"_",1:n)
#   object@reductions[[reduction_name]]@cell.embeddings <- embeddings
#   object@reductions[[reduction_name]]@key <- ifelse(test = grepl("_$",key_name), yes = key_name, no = paste0(key_name,"_"))
#   return(object)
# }

#------- process seurat to loom


# library(loomR)
# library(hdf5r)
# library(Seurat)
# setwd('/Users/jsjoyal/Desktop/SCAP/test_data/projects/')
# 
# seur <- readRDS("/Users/jsjoyal/Desktop/SCAP/test_data/GSE129516.shortname.obj.rds")
# new_dir <- 'GSE129516_shortname'
# 
# dir.create(new_dir)
# project_dir <- paste0(new_dir,'/')
# 
# down_sample <- round(6000/length(unique(seur$seurat_clusters)))
# seur <- subset(seur, downsample = down_sample, idents = 'seurat_clusters')
# # if(dim(seur)[2]>6000){
# #   seur <- subset(seur, cells = sample(Cells(seur), 6000))
# # }
# 
# assays <- names(seur@assays)
# 
# flag <- FALSE
# if(any(assays == 'integrated')){
#   assays <- assays[-grep('integrated|SCT',assays)]
# }else if(any(assays=='SCT')){
#     assays <- assays[-grep('integrated|RNA',assays)]
#     flag <- TRUE
# }
# 
# for(assay in assays){
#   DefaultAssay(seur) <- assay
#   if(flag == TRUE && assay == "SCT"){
#     pseudo.assay <- "RNA"
#   }else{
#     pseudo.assay <- assay
#   }
#   filename <- paste0(project_dir,pseudo.assay,".loom")
#   loomR::create(filename = filename, data = seur[[assay]]@data, calc.count = F, overwrite = T)
#   data <- loomR::connect(filename = filename, mode = "r+")
#   data$link_delete('row_attrs/Gene')
#   
#   # add metadata
#   meta.data <- seur@meta.data
#   colnames(meta.data) <- paste0(colnames(meta.data),"_meta_data")
# 
#   data$add.col.attribute(as.list(meta.data))
#   data$add.row.attribute(list(features = rownames(seur[[pseudo.assay]])))
# 
#   # add reduction embeddings
#   reduction_names <- names(seur@reductions)
#   for(i in 1:length(reduction_names)){
#     reductions <- as.data.frame(seur@reductions[[i]]@cell.embeddings)
#     assay_used <- tolower(seur@reductions[[i]]@assay.used)
#     if(!grepl(paste0('(?![a-z])(?<![a-z])(',assay_used,')'),reduction_names[i],perl=T,ignore.case = T)){
#       reduction_names[i] <- paste0(reduction_names[i],'_',assay_used)
#     }
#     n <- ncol(reductions)
#     if(n>3){
#       for(j in 2:3){
#         tmp <- reductions[,1:j]
#         tmp_name <- reduction_names[i]
#         if(!grepl(paste0('(?![a-z])(?<![a-z])(',j,'d)'),tmp_name,perl=T,ignore.case = T)){
#           tmp_name <- paste0(tmp_name,'_',j,'d')
#         }
#         if(grepl(' ', tmp_name)){
#           tmp_name <- gsub(' ','_',tmp_name)
#         }
#         colnames(tmp) <- paste0(tmp_name,'_',1:j,'_reduction')
#         data$add.col.attribute(as.list(tmp))
#       }
#     }else{
#       if(!grepl(paste0('(?![a-z])(?<![a-z])(',n,'d)'),reduction_names[i],perl=T,ignore.case = T)){
#         reduction_names[i] <- paste0(reduction_names[i],'_',n,'d')
#       }
#       if(grepl(' ', reduction_names[i])){
#         reduction_names[i] <- gsub(' ','_',reduction_names[i])
#       }
#       colnames(reductions) <- paste0(reduction_names[i],'_',1:n,'_reduction')
#       data$add.col.attribute(as.list(reductions))
#     }
#   }
#   data$close_all()
# }
# 
