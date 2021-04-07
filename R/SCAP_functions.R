library("shiny")
library("cowplot")
library("ggplot2")
library("ggridges")
library("ggthemes")
library("gtools")
library("loomR")
library("Matrix")
library("plotly")
library("Seurat")
library("viridis")
library("Nebulosa")
library("Seurat")
library("SeuratDisk")
library("reticulate")

palette <- colorRampPalette(c("lightgrey", viridis(10)))

check_if_obs_cat <- function(obs_df){
  len <- nrow(obs_df)
  uni <- apply(obs_df, 2, function(x){length(unique(x))})
  cat <- (uni > (0.05*len))
  num_check <- apply(obs_df, 2, function(x){suppressWarnings(all(!is.na(as.numeric(as.character(x)))))})
  return(cat*num_check == 0) # TRUE = discrete, FALSE = continuous
}

save_figure <- function(file_type, file_name, units, height, width, resolution, to_plot){
  if(file_type == "png"){
    png(filename = file_name, units = units, height = height, width = width, res = resolution)
    to_plot
    dev.off()
  }else if(file_type == "jpeg"){
    jpeg(filename = file_name, units = units, height = height, width = width, res = resolution)
    to_plot
    dev.off()
  }else if(file_type == "tiff"){
    tiff(filename = file_name, units = units, height = height, width = width, res = resolution)
    to_plot
    dev.off()
  }else if(file_type == "pdf"){
    pdf(filename = file_name, units = units, height = height, width = width, res = resolution)
    to_plot
    dev.off()
  }
}

loomToSeurat <- function(obj, loom, dir, file){
  thisDate <- gsub(":","_",gsub(" ","_",date()))
  print(paste("thisDate:", thisDate))
  print(paste("obj:", obj))
  print(paste("dir:", dir))
  print(paste("file:", file))
  print("loom:")
  print(str(loom))
  seur <- readRDS(obj)
  print(nrow(seur@meta.data))
  print("seur meta data:")
  print(str(seur@meta.data))
  print('loom # cells:')
  print(loom[[1]]$col.attrs[[1]]$dims)
  if(nrow(seur@meta.data) != loom[[1]]$col.attrs[[1]]$dims){
    showModal(modalDialog(p("Error: Seurat object and loom object have different number of cells and likely do not coresspond to the same experiment"), title = paste0("ERROR")), session = getDefaultReactiveDomain())
    return(-1)
  }
  loom_metaData <- loom[[1]]$get.attribute.df(attributes = names(loom[[1]]$col.attrs)[grep("_meta_data$", names(loom[[1]]$col.attrs))])
  print("loom_metaData:")
  print(colnames(loom_metaData))
  colnames(loom_metaData) <- sub("_meta_data$","",colnames(loom_metaData))
  print("loom_metaData:")
  print(colnames(loom_metaData))
  same_names <- colnames(loom_metaData)[which(colnames(loom_metaData)%in%colnames(seur@meta.data))]
  if(length(same_names)>0){
    for(i in 1:length(same_names)){
      if(identical(loom_metaData[,same_names[i]],seur@meta.data[,same_names[i]])){
        loom_metaData <- loom_metaData[,-which(colnames(loom_metaData)==same_names[i])]
      }else{
        colnames(loom_metaData)[which(colnames(loom_metaData)==same_names[i])] <- paste0(same_names[i],"_",thisDate)
      }
    }
  }
  print(paste("loom_metaData:", colnames(loom_metaData)))
  seur@meta.data <- cbind(seur@meta.data, loom_metaData)
  print(paste('new mata data:'))
  print(str(seur@meta.data))
  tryCatch(
    {
      seur <- UpdateSeuratObject(seur)
    },
    error=function(cond){
      message('Error: Could not update Seurat Object:')
      message(cond)
    },
    warning=function(cond){
      message('Warning when trying to update Seurat Object:')
      message(cond)
      seur <- UpdateSeuratObject(seur)
    },
    finally={
      message('Saving Seurat Object...')
    }
  )
  if(is.null(file)){
    saveRDS(seur, obj)
  }else{
    saveRDS(seur, paste0(dir,"/",file))
  }
  return(1)
}

loomToSeurat_legacy <- function(obj, loom){
  thisDate <- gsub(":","_",gsub(" ","_",date()))
  # print(paste("thisDate:", thisDate))
  # print(paste("obj:", obj))
  # print(paste("dir:", dir))
  # print(paste("file:", file))
  # print("loom:")
  # print(str(loom))
  print(nrow(obj@meta.data))
  print("obj meta data:")
  print(str(obj@meta.data))
  print('loom # cells:')
  print(loom[[1]]$col.attrs[[1]]$dims)
  if(nrow(obj@meta.data) != loom[[1]]$col.attrs[[1]]$dims){
    showModal(modalDialog(p("Error: Seurat object and loom object have different number of cells and likely do not coresspond to the same experiment"), title = paste0("ERROR")), session = getDefaultReactiveDomain())
    return(-1)
  }
  loom_metaData <- loom[[1]]$get.attribute.df(attributes = names(loom[[1]]$col.attrs)[grep("_meta_data$", names(loom[[1]]$col.attrs))])
  print("loom_metaData:")
  print(colnames(loom_metaData))
  colnames(loom_metaData) <- sub("_meta_data$","",colnames(loom_metaData))
  print("loom_metaData:")
  print(colnames(loom_metaData))
  same_names <- colnames(loom_metaData)[which(colnames(loom_metaData)%in%colnames(obj@meta.data))]
  if(length(same_names)>0){
    for(i in 1:length(same_names)){
      if(identical(loom_metaData[,same_names[i]],obj@meta.data[,same_names[i]])){
        loom_metaData <- loom_metaData[,-which(colnames(loom_metaData)==same_names[i])]
      }else{
        colnames(loom_metaData)[which(colnames(loom_metaData)==same_names[i])] <- paste0(same_names[i],"_",thisDate)
      }
    }
  }
  print(paste("loom_metaData:", colnames(loom_metaData)))
  obj@meta.data <- cbind(obj@meta.data, loom_metaData)
  print(paste('new mata data:'))
  print(str(obj@meta.data))
  # if(is.null(file)){
  #   saveRDS(seur, obj)
  # }else{
  #   saveRDS(seur, paste0(dir,"/",file))
  # }
  return(obj)
}


seuratToLoom <- function(obj, dir){
  #library(loomR)
  #library(hdf5r)
  #library(Seurat)

  ## Check file permissions before attempting to read the file
  file_access <- file.access(obj, mode = 4)
  if(file_access == -1){
	  showNotification('Error: Permission denied! Make sure you have changed file permissions for this object!', type = 'error')
	  return(0) ## Return failed error code
  }

  ## Make sure the object is really .rds before attempting to read it
  if(endsWith(obj,".rds")){
    seur <- try(readRDS(obj))
  }else if(grepl(".h5|.h5ad",obj) ){
    seur <- try(ReadH5AD(obj)) ## Read in scanpy object and convert to
    if(grepl('^seurat',class(seur)[1],ignore.case = T)){
      Idents(seur) <- "louvain_scanpy_clusters"
      seur@assays$RNA@counts <- seur@assays$RNA@data
    }
  }

  if(class(seur) == "try-error"){
    showNotification('Error: The selected file is labeled .rds but is not actually a valid .rds file!', type = 'error')
    return(0) ## Return failed error code
  }
  
  if(!grepl('^seurat',class(seur)[1],ignore.case = T)){
    showNotification('Error: The selected object is of class ', class(seur)[1], ' but must be of class Seurat', type = 'error')
  }
  
  current_version <- 3
  if(strsplit(as.character(seur@version), split = '\\.')[[1]][1]<current_version){
    showNotification(paste0("Warning: The selected seurat object is out of date (version: ", seur@version, "). Updating to object now..."), type = 'error')
    seur <- UpdateSeuratObject(seur)
    showNotification(paste0("Update complete"), type = 'message')
  }
  
  # make sure percent.mito is added to the object
  if(any(names(seur@assays)=='RNA')){
    DefaultAssay(seur) <- 'RNA'
  }
  if(any(colnames(seur@meta.data) == 'percent.mito') == FALSE){
    if(any(grepl('^MT-', rownames(seur)))){
      seur$percent.mito <- PercentageFeatureSet(seur, pattern = '^MT-')
    }else if(any(grepl('^mt-', rownames(seur)))){
      seur$percent.mito <- PercentageFeatureSet(seur, pattern = '^mt-')
    }else if(any(grepl("^Mt-", rownames(seur)))){
      seur$percent.mito <- PercentageFeatureSet(seur, pattern = '^Mt-')
    }else{
      seur$percent.mito <- 0
      showNotification(paste0("Mitochondrial genes were not identified."), type = 'warning')
    }
  }
  
  project_dir <- paste0(dir,'/')
  
  #down_sample <- round(6000/length(unique(seur$seurat_clusters)))
  #seur <- subset(seur, downsample = down_sample, idents = 'seurat_clusters')
  # if(dim(seur)[2]>6000){
  #   seur <- subset(seur, cells = sample(Cells(seur), 6000))
  # }
  
  assays <- names(seur@assays)
  
  if(length(assays)>1){
    dims <- lapply(seur@assays, function(x){
      return(dim(x@data))
    })
    names(dims) <- assays
    assay_index <- NULL
    if(any(names(assays)=="RNA" & dims[["RNA"]][2] > 0)){
      for(i in 1:length(dims)){
        if(dims[[i]][2]!=dims[["RNA"]][2]){
          showNotification(paste0("Warning: the ", assays[i], "assay does not have the same number of cells (",dims[[i]][2],") as the RNA assay (", dims[["RNA"]][2],") and will be removed."), type = 'warning')
          assay_index <- c(assay_index,i)
        }
      }
    }else{
      first_assay <- dims[[1]][2]
      for(i in 2:length(dims)){
        if(first_assay!=dims[[i]][2]){
          showNotification(paste0("Warning: the ", assays[i], "assay does not have the same number of cells (",dims[[i]][2],") as the reference assay (", assay[1],": ",dims[["RNA"]][2],") and will be removed."), type = 'warning')
          assay_index <- c(assay_index,i)
        }
      }
    }
    if(!is.null(assay_index)){
      assays <- assays[-assay_index] 
    }
  }
  
  for(assay in assays){
    DefaultAssay(seur) <- assay
    
    filename <- paste0(project_dir,assay,".loom")
    loomR::create(filename = filename, data = seur[[assay]]@data, calc.count = F, overwrite = T)
    data <- loomR::connect(filename = filename, mode = "r+")
    data$link_delete('row_attrs/Gene')
    
    # add metadata
    meta.data <- seur@meta.data
    colnames(meta.data) <- paste0(colnames(meta.data),"_meta_data")
    
    data$add.col.attribute(as.list(meta.data))
    data$add.row.attribute(list(features = rownames(seur[[assay]])))
    
    # add reduction embeddings
    reduction_names <- names(seur@reductions)
    for(i in 1:length(reduction_names)){
      reductions <- as.data.frame(seur@reductions[[i]]@cell.embeddings)
      if(nrow(reductions)==0){
        showNotification(paste0("Warning: There was no data found for the ", reduction_names[i], " reduction. Skipping this reduction."), type = 'warning')
        next
      }
      assay_used <- tolower(seur@reductions[[i]]@assay.used)
      #if(!grepl(paste0('(?![a-z])(?<![a-z])(',assay_used,')'),reduction_names[i],perl = T,ignore.case = T)){
      if(!grepl(assay_used,reduction_names[i], ignore.case = T)){
        reduction_names[i] <- paste0(reduction_names[i],'_',assay_used)
      }
      n <- ncol(reductions)
      if(n>3){
        for(j in 2:3){
          tmp <- reductions[,1:j]
          tmp_name <- reduction_names[i]
          if(!grepl(paste0('(?![a-z])(?<![a-z])(',j,'d)'),tmp_name,perl=T,ignore.case = T)){
            tmp_name <- paste0(tmp_name,'_',j,'d')
          }
          if(grepl(' ', tmp_name)){
            tmp_name <- gsub(' ','_',tmp_name)
          }
          colnames(tmp) <- paste0(tmp_name,'_',1:j,'_reduction')
          data$add.col.attribute(as.list(tmp))
        }
      }else{
        if(!grepl(paste0('(?![a-z])(?<![a-z])(',n,'d)'),reduction_names[i],perl=T,ignore.case = T)){
          reduction_names[i] <- paste0(reduction_names[i],'_',n,'d')
        }
        if(grepl(' ', reduction_names[i])){
          reduction_names[i] <- gsub(' ','_',reduction_names[i])
        }
        colnames(reductions) <- paste0(reduction_names[i],'_',1:n,'_reduction')
        data$add.col.attribute(as.list(reductions))
      }
    }
    data$close_all()
  }
  return(1)
}

try_seurat_update <- function(seur){
  tryCatch(
    {
      seur <- UpdateSeuratObject(seur)
    },
    error=function(cond){
      message('Error: Could not update Seurat Object:')
      message(cond)
    },
    warning=function(cond){
      message('Warning when trying to update Seurat Object:')
      message(cond)
      seur <- UpdateSeuratObject(seur)
    },
    finally={
      return(seur)
    }
  )
}

detect_legacy_scap <- function(loom){
  col <- names(loom[[1]][['col_attrs']])
  if(any(grepl("_meta_data$", col)) && any(grepl("_reduction$", col))){
    return(1)
  }else{
    return(0)
  }
}

scap_to_h5ad <- function(in_file, out_path, old_file, modality = "RNA"){
  anndata <- import('anndata')
  if(grepl("\\.rds$", in_file, ignore.case = TRUE)){
    obj <- readRDS(in_file)
  }else if(grepl("\\.loom$", in_file, ignore.case = TRUE)){
    obj <- connect(in_file)
  }else{
    return(0)
  }
  if(identical(class(obj)[1], "Seurat")){
    obj <- UpdateSeuratObject(obj)
    message("Converting to H5Seurat")
    SaveH5Seurat(obj, filename = sub("\\.h5ad$", ".h5Seurat", out_path))
    message("Converting to H5ad")
    Convert(sub("\\.h5ad$", ".h5Seurat", out_path), dest = "h5ad", assay = modality)
  }else if(identical(class(obj)[1], "loom")){
    legacy <- detect_legacy_scap(list(obj))
    if(legacy){
      if(is.null(old_file)){
        return(-1)
      }else{
        old_obj <- readRDS(old_file)
        if(identical(class(old_obj)[1], 'Seurat')){
          seur <- loomToSeurat_legacy(old_obj, list(obj))
          if(identical(seur, -1)){
            return(-2)
          }else{
            seur <- UpdateSeuratObject(seur)
            message("Converting to H5Seurat")
            SaveH5Seurat(seur, filename = sub("\\.h5ad$", ".h5Seurat", out_path))
            message("Converting to H5ad")
            Convert(sub("\\.h5ad$", ".h5Seurat", out_path), dest = "h5ad", assay = modality)
          }
        }else{
          return(-3)
        }
      }
    }else{
      lf <- anndata$read_loom(file)
      lf$write(out_path)
    }
  }else{
    return(-4)
  }
  #------ required to fix anndata.write ValueError: '_index' is a reserved name for dataframe columns. Above error raised while writing key 'raw/var' of <class 'h5py._hl.files.File'> from /. ------#
  cat(file = stderr(), "raw/var correction...\n")
  a = anndata$read(out_path)
  py$tmp = a
  py_run_string("tmp.__dict__['_raw'].__dict__['_var'] = tmp.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})")
  a = py$tmp
  a$write(out_path)
  #--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
  return(1)
}



reorder_levels <- function(x){
  if(!is.factor(x)){
    x <- as.factor(x)
  }
  factor(x,mixedsort(levels(x)))
}

seurat_reduc_key_name <- function(object, reduction_name = NULL, key_as_name = TRUE, key_name = NULL){
  if(reduction_name==FALSE) stop('ERROR: you must specify the reduction name')
  if(key_as_name == FALSE && is.null(key_name)){
    stop('ERROR: If you do not want to set the key name to the reduction name, a key name must be specified')
  }else if(key_as_name == TRUE && is.null(key_name)){
    print('setting key_name to reduction_name...')
    key_name <- reduction_name
  }
  embeddings <- object@reductions[[reduction_name]]@cell.embeddings
  n <- ncol(embeddings)
  colnames(embeddings) <- paste0(key_name,"_",1:n)
  object@reductions[[reduction_name]]@cell.embeddings <- embeddings
  object@reductions[[reduction_name]]@key <- ifelse(test = grepl("_$",key_name), yes = key_name, no = paste0(key_name,"_"))
  return(object)
}

######### Programatic trigger ###########
makeReactiveTrigger <- function() {
  rv <- reactiveValues(a = 0)
  list(
    depend = function() {
      rv$a
      invisible()
    },
    trigger = function() {
      rv$a <- isolate(rv$a + 1)
    }
  )
}
###### reduc_key ##########
reduc_key <- function(key){
  if(!any(key == c("TSNE", "UMAP", "PCA", "DIFF", "CCA"))){
    showNotification("Reduction unknown. Axis labels may be poorly formatted")
    return(key)
  }
  if(key == "PCA")
    return("PC")
  if(key == "DIFF")
    return("DC")
  if(key == "TSNE")
    return("tSNE")
  if(key == "UMAP")
    return("UMAP")
  if(key == "CCA")
    return("CCA")
}

###### percentAbove #######
percentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}

######### dotPlot ############
dotPlot <- function(
  data.features = NULL,
  assay = NULL,
  features = NULL,
  cols = c("blue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = 'radius',
  scale.min = NA,
  scale.max = NA,
  ...
) {
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  
  if(is.null(features)) return(NULL)

  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = percentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      data.use <- scale(x = data.use)
      data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = rev(x = features)
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  color.by <- ifelse(test = is.null(x = split.by), yes = 'avg.exp.scaled', no = 'colors')
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }
  
  ax.x <- list(
    title = "Features",
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = TRUE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='outside'
  )
  ax.y <- list(
    title = "Identity",
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = TRUE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='outside'
  )
  
  data.plot$id <- reorder_levels(data.plot$id)
  plot <- ggplot(data = data.plot, mapping = aes(x = features.plot, y = id)) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
    #scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(
      x = 'Features',
      y = 'Identity'
    )
  if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }

  plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
  
  return(plot)
}

######## dimPlotlyOutput #######
dimPlotlyOutput <- function(assay.in, reduc.in, group.by, annot_panel = NULL, tmp_annotations = NULL, low.res){

  n <- length(colnames(reduc.in))
  if(n<2 || n>3){
    showNotification(ui = paste0('Error: Invalid number of dimensions for selected reduction: ',n),type = 'error')
  }

  group.by[[1]] <- reorder_levels(group.by[[1]])

  plot.data <- cbind(reduc.in, group.by[[1]])
  colnames(plot.data) <- c(colnames(reduc.in), names(group.by)[1])
  rownames(plot.data) <- names(group.by[[1]])

  ax.x <- list(
    title = colnames(reduc.in)[1],
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'#,
    #range = c(xlimits[1],xlimits[2])
  )
  ax.y <- list(
    title = colnames(reduc.in)[2],
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'#,
    #range = c(ylimits[1],ylimits[2])
  )
  if(length(colnames(reduc.in)) == 3){
    ax.z <- list(
      title = colnames(reduc.in)[3],
      zeroline = FALSE,
      showline = TRUE,
      showticklabels = FALSE,
      showgrid = FALSE,
      mirror=TRUE,
      ticks='none'#,
      #range = c(zlimits[1],zlimits[2])
    )
  }
  
  col <- if(annot_panel == 'cell_annotation_custom'){
    tmp_annotations
  }else if(n==2){
    plot.data[,3]
  }else{
    plot.data[,4]
  }

  if(n == 2){
    p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2], 
                 type = 'scatter', mode = 'markers', key = ~rownames(plot.data), alpha = 0.6, stroke = I('dimgrey'),
                 color = col, text =  ~paste0(
                   names(group.by)[1],": ", plot.data[,3],"\n</br>",
                   colnames(reduc.in)[1],": ", format(plot.data[,1],digits=3),"\n",
                   "</br>",colnames(reduc.in)[2],": ", format(plot.data[,2],digits=3)),
                 hovertemplate = paste0('<b>%{text}</b><extra></extra>')
    ) %>% layout(title = ifelse(test = annot_panel == 'cell_annotation_custom', yes = paste0(assay.in, " data coloured by custom annotations"), no = paste0(assay.in, " data coloured by ", names(group.by))) ,xaxis = ax.x, yaxis = ax.y, dragmode='lasso')
  }else{
    p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2], z = plot.data[,3],
                 type = 'scatter3d', mode = 'markers', key = ~rownames(plot.data), stroke = I('dimgrey'),
                 color = col, text =  ~paste0(
                   names(group.by)[1],": ", plot.data[,4],"\n</br>",
                   colnames(reduc.in)[1],": ", format(plot.data[,1],digits=3),"\n",
                   "</br>",colnames(reduc.in)[2],": ", format(plot.data[,2],digits=3), "\n",
                   "</br>",colnames(reduc.in)[3],": ", format(plot.data[,3],digits=3)), 
                 hovertemplate = paste0('<b>%{text}</b><extra></extra>')
    ) %>% layout(title = ifelse(test = annot_panel == 'cell_annotation_custom', yes = paste0(assay.in, " data coloured by custom annotations"), no = paste0(assay.in, " data coloured by ", names(group.by))) ,scene = list(xaxis = ax.x, yaxis = ax.y, zaxis = ax.z, dragmode='lasso'),legend = list(x = 100, y = 0.5))
  }

  if(low.res == 'yes'){
    return(p %>% toWebGL())
  }else{
    return(p)
  }
}

####### featurePlotlyOutput ##########
featurePlotlyOutput <- function(assay.in, reduc.in, group.by, feature.in, low.res){
  
  n <- length(colnames(reduc.in))
  if(n<2 || n>3){
    showNotification(ui = paste0('Error: Invalid number of dimensions for selected reduction: ',n),type = 'error')
  }

  if(ncol(feature.in)==0) return(NULL)
    
  group.by[[1]] <- reorder_levels(group.by[[1]])

  plot.data <- do.call(cbind, list(reduc.in, group.by, feature.in))
  rownames(plot.data) <- names(group.by[[1]])

  ax.x <- list(
    title = colnames(reduc.in)[1],
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'#,
    #range = c(xlimits[1],xlimits[2])
  )
  ax.y <- list(
    title = colnames(reduc.in)[2],
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'#,
    #range = c(ylimits[1],ylimits[2])
  )
  if(length(colnames(reduc.in)) == 3){
    ax.z <- list(
      title = colnames(reduc.in)[3],
      zeroline = FALSE,
      showline = TRUE,
      showticklabels = FALSE,
      showgrid = FALSE,
      mirror=TRUE,
      ticks='none'#,
      #range = c(zlimits[1],zlimits[2])
    )
  } 

  if(n == 2){
    p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2],
                 type = 'scatter', mode = 'markers', key = ~rownames(plot.data),
                 color = plot.data[,4], colors = palette(100), opacity = 0.6, text =  ~paste0(
                   names(group.by)[1],": ", plot.data[,3],"\n",
                   "</br>",'Expression',": ", format(plot.data[,4],digits=3),"\n",
                   "</br>",colnames(reduc.in)[1],"_1: ", format(plot.data[,1],digits=3),"\n",
                   "</br>",colnames(reduc.in)[2],"_2: ", format(plot.data[,2],digits=3), "\n"), 
                 marker = list(size = 5),
                 hovertemplate = paste0('<b>%{text}</b>',
                                        '<extra></extra>')
    ) %>% layout(title = feature.in ,xaxis = ax.x, yaxis = ax.y, dragmode='lasso')
  }else{
    p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2], z = plot.data[,3],
                 type = 'scatter3d', mode = 'markers', key = ~rownames(plot.data),
                 color = plot.data[,5], colors = palette(100), opacity = 0.6, text =  ~paste0(
                   names(group.by)[1],": ", plot.data[,4],"\n",
                   "</br>",'Expression',": ", format(plot.data[,5],digits=3),"\n",
                   "</br>",colnames(reduc.in)[1],"_1: ", format(plot.data[,1],digits=3),"\n",
                   "</br>",colnames(reduc.in)[2],"_2: ", format(plot.data[,2],digits=3), "\n",
                   "</br>",colnames(reduc.in)[3],"_3: ", format(plot.data[,3],digits=3), "\n"), 
                 marker = list(size = 3),
                 hovertemplate = paste0('<b>%{text}</b>',
                                        '<extra></extra>')
    ) %>% layout(title = feature.in , scene = list(xaxis = ax.x, yaxis = ax.y, zaxis = ax.z), dragmode='lasso' ,legend = list(x = 100, y = 0.5))
  }
  #print(paste("featureplot plot:", Sys.time()-t4))
  #print(paste("featureplot total:", Sys.time()-t1))
  if(low.res == 'yes'){
    return(p%>% toWebGL())
  }else{
    return(p)
  }
}

#### 
## Nebulosa version
featurePlotlyOutput_nebulosa <-  function(assay.in, reduc.in, group.by, feature.in, low.res){
  if(ncol(feature.in) < 1) return(NULL)
  data.features = as.matrix(t(feature.in))
  data_seurat <- CreateSeuratObject(counts = data.features,
                                    meta.data = group.by,
                                    min.cells = 0,
                                    min.features = 0)
  
  data_seurat[['umap']] <- CreateDimReducObject(embeddings = as.matrix(reduc.in),
                                    assay = "RNA",
                                    key = unique(sub("_[0-9]$", "", colnames(reduc.in)))) 
  
  if(ncol(feature.in) == 1){
    return(plot_density(data_seurat, features = colnames(feature.in), size = 2) + theme_cowplot())
  }else{
    p <- plot_density(data_seurat, features = colnames(feature.in), joint = TRUE, combine = FALSE, size = 2) 
    return(p[[length(p)]] + theme_cowplot())
  }
}

####

split_dot_plot <- function(data.features = NULL, 
                           assay = NULL, 
                           features = NULL, 
                           split.by = NULL,
                           col.max = 2.5,
                           col.min = -2.5,
                           cols = c("blue","red"),
                           dot.min = 0,
                           dot.scale = 6,
                           scale.by = "radius",
                           scale.max = NA,
                           scale.min = NA){
  
  if(is.null(features)) return(NULL)
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- split.by[[1]]
    cols <- cols[1:length(x = unique(x = splits))]
    names(x = cols) <- unique(x = splits)
    data.features$id <- paste(data.features$id, splits, sep = 'SPLITBYTAG')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "SPLITBYTAG", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = percentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      data.use <- scale(x = data.use)
      data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (!is.null(x = split.by)) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = rev(x = features)
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  
  data.plot$split <- sub(".*SPLITBYTAG","",data.plot$id)
  data.plot$id <- sub("SPLITBYTAG.*","",data.plot$id)
  
  data.plot <- data.plot
  data.plot$id <- as.factor(data.plot$id)
  data.plot$split <- as.factor(data.plot$split)
  
  data.plot$id <- reorder_levels(data.plot$id)
  
  labels <- data.frame(features.plot = rep(levels(data.plot$features.plot),length(levels(data.plot$split))))
  labels$id <- rep(levels(data.plot$id)[length(levels(data.plot$id))], nrow(labels))
  labels <- labels[order(labels$features.plot),]
  labels$split <- rep(levels(data.plot$split),nrow(labels)/length(levels(data.plot$split)))
  
  labels$features.plot <- as.factor(labels$features.plot)
  labels$id <- as.factor(labels$id)
  labels$split <- as.factor(labels$split)
  labels$y <- as.factor('ZZZZZZ') #length(levels(data.plot$id))+1
  
  p <- ggplot(data = data.plot, mapping = aes(x = features.plot, y = id)) + 
    geom_point(aes(color = avg.exp.scaled, size = pct.exp, group = split), position = position_dodge(1))  +
    geom_text(data = labels, aes(x = features.plot, y = y, label = split, group = split), position = position_dodge(width = 1), size=5) + 
    scale_y_discrete(breaks = as.factor(levels(data.plot$id))) +
    scale_color_gradient(low = 'blue', high = 'red') + 
    xlab('Features') +
    ylab('Identity') +
    theme_base()
  
  #print('layer_1:')
  #print(layer_data(p,1))
  layer <- layer_data(p,2)
  layer <- layer[order(layer$x),]
  layer$y <- layer$y-0.5
  #print('layer:')
  #print(layer)
  n_groups <- length(unique(layer$group))
  start <- 1
  for(i in 1:nrow(layer)){
    if(start%%n_groups==1){
      layer$x[i] <- layer$x[i] - 0.15
    }
    if(start%%n_groups==0){
      layer$x[i] <- layer$x[i] + 0.15
    }
    if(start==n_groups){
      start <- 1
    }else{
      start <- start+1
    }
  }
  
  cnt <- 1
  #print('l_data:')
  while(cnt<1000){
    l_data <- layer[cnt:(cnt+n_groups-1),]
    #print(l_data)
    p <- p + geom_line(data = l_data, aes(x = x, y = y))
    cnt <- cnt + n_groups
    if(cnt>(length(levels(data.plot$features.plot))*n_groups)){
      break
    }
  }
  return(p)
}

FacetTheme <- function(...) {
  return(theme(
    strip.background = element_blank(),
    strip.text = element_text(face = 'bold'),
    # Validate the theme
    validate = TRUE,
    ...
  ))
}

ggRidgePlot <- function(
  data,
  sort = FALSE,
  same.y.lims = TRUE,
  adjust = 1,
  pt.size = 0,
  cols = NULL,
  log = FALSE,
  fill.by = NULL,
  flip = NULL
) {

  cat(file = stderr(), "data$expression: ", summary(data$expression), "\n")

  if (!is.data.frame(x = data) || ncol(x = data) < 2) {
    stop("RidgePlotly requires a data frame with >1 column")
  }

  if (sort) {
    data$feature <- as.vector(x = data$feature)
    data$ident <- as.vector(x = data$ident)
    # build matrix of average expression (#-features by #-idents), lexical ordering
    avgs.matrix <- sapply(
      X = split(x = data, f = data$ident),
      FUN = function(df) {
        return(tapply(
          X = df$expression,
          INDEX = df$feature,
          FUN = mean
        ))
      }
    )
    idents.order <- hclust(d = dist(x = t(x = L2Norm(mat = avgs.matrix, MARGIN = 2))))$order
    avgs.matrix <- avgs.matrix[,idents.order]
    avgs.matrix <- L2Norm(mat = avgs.matrix, MARGIN = 1)
    # order feature clusters by position of their "rank-1 idents"
    position <- apply(X = avgs.matrix, MARGIN = 1, FUN = which.max)
    mat <- hclust(d = dist(x = avgs.matrix))$merge
    orderings <- list()
    for (i in 1:nrow(mat)) {
      x <- if (mat[i,1] < 0) -mat[i,1] else orderings[[mat[i,1]]]
      y <- if (mat[i,2] < 0) -mat[i,2] else orderings[[mat[i,2]]]
      x.pos <- min(x = position[x])
      y.pos <- min(x = position[y])
      orderings[[i]] <- if (x.pos < y.pos) {
        c(x, y)
      } else {
        c(y, x)
      }
    }
    features.order <- orderings[[length(x = orderings)]]
    data$feature <- factor(
      x = data$feature,
      levels = unique(x = sort(x = data$feature))[features.order]
    )
    data$ident <- factor(
      x = data$ident,
      levels = unique(x = sort(x = data$ident))[rev(x = idents.order)]
    )
  } else {
    data$feature <- factor(x = data$feature, levels = unique(x = data$feature))
  }
  # if (log) {
  #   noise <- rnorm(n = nrow(x = data)) / 200
  #   data$expression <- data$expression + 1
  # } else {
  #   noise <- rnorm(n = nrow(x = data)) / 100000
  # }
  # for (f in unique(x = data$feature)) {
  #   if (all(data$expression[(data$feature == f)] == data$expression[(data$feature == f)][1])) {
  #     warning(
  #       "All cells have the same value of ",
  #       f,
  #       call. = FALSE,
  #       immediate. = TRUE
  #     )
  #   } else {
  #     data$expression[(data$feature == f)] <- data$expression[(data$feature == f)] + noise[(data$feature == f)]
  #   }
  # }

  # if (flip) {
  #   x <- 'ident'
  #   x.label <- 'Identity'
  #   y <- 'expression'
  #   y.label <- 'Expression Level'
  # } else {
  #   y <- 'ident'
  #   y.label <- 'Identity'
  #   x <- 'expression'
  #   x.label <- 'Expression Level'
  # }

  data$ident <- reorder_levels(data$ident)
  
  plot <- ggplot(data, aes(y=ident, x=expression, fill=ident)) +
  geom_density_ridges(alpha=0.6, bandwidth=0.5) +
  xlab("Expression") +
  ylab("Identity") + 
  facet_grid(. ~ feature, scales = (if (same.y.lims) 'fixed' else 'free')) +
  FacetTheme(
    panel.spacing = unit(0, 'lines'),
    panel.background = element_rect(fill = NA, color = "black"),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y.right = element_text(angle = 0)
  )

  return(plot)
}

