#' Functions For Shiny App 
#' 
#' @import cowplot
#' @import dplyr
#' @import ggplot2
#' @import ggthemes
#' @import gtools
#' @import hdf5r
#' @import loomR
#' @import Matrix
#' @import MODIS
#' @import plotly
#' @import presto
#' @import Seurat
#' @import shiny
#' @import shinycssloaders
#' @import shinyFiles
#' @import shinyjqui
#' @import shinythemes

#my.pal <- c(RColorBrewer::brewer.pal(12,'Set3'),RColorBrewer::brewer.pal(12,'Paired'))

seuratToLoom <- function(obj, dir){
  #library(loomR)
  #library(hdf5r)
  #library(Seurat)
  
  seur <- readRDS(obj)
  
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
    }
  }
  
  project_dir <- paste0(dir,'/')
  
  #down_sample <- round(6000/length(unique(seur$seurat_clusters)))
  #seur <- subset(seur, downsample = down_sample, idents = 'seurat_clusters')
  # if(dim(seur)[2]>6000){
  #   seur <- subset(seur, cells = sample(Cells(seur), 6000))
  # }
  
  assays <- names(seur@assays)
  
  flag <- FALSE
  if(any(assays == 'integrated')){
    assays <- assays[-grep('integrated|SCT',assays)]
  }else if(any(assays=='SCT')){
    assays <- assays[-grep('integrated|RNA',assays)]
    flag <- TRUE
  }
  
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
    
    # add metadata
    meta.data <- seur@meta.data
    colnames(meta.data) <- paste0(colnames(meta.data),"_meta_data")
    
    data$add.col.attribute(as.list(meta.data))
    data$add.row.attribute(list(features = rownames(seur[[pseudo.assay]])))
    
    # add reduction embeddings
    reduction_names <- names(seur@reductions)
    for(i in 1:length(reduction_names)){
      reductions <- as.data.frame(seur@reductions[[i]]@cell.embeddings)
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

withBusyIndicatorCSS <- "
.btn-loading-container {
margin-left: 10px;
font-size: 1.2em;
}
.btn-done-indicator {
color: green;
}
.btn-err {
margin-top: 10px;
color: red;
}
"

withBusyIndicatorUI <- function(button) {
  id <- button[['attribs']][['id']]
  div(
    shinyjs::useShinyjs(),
    singleton(tags$head(
      tags$style(withBusyIndicatorCSS)
    )),
    `data-for-btn` = id,
    button,
    span(
      class = "btn-loading-container",
      shinyjs::hidden(
        icon("spinner", class = "btn-loading-indicator fa-spin"),
        icon("check", class = "btn-done-indicator")
      )
    ),
    shinyjs::hidden(
      div(class = "btn-err",
          div(icon("exclamation-circle"),
              tags$b("Error: "),
              span(class = "btn-err-msg")
          )
      )
    )
  )
}
withBusyIndicatorServer <- function(buttonId, expr) {
  # UX stuff: show the "busy" message, hide the other messages, disable the button
  loadingEl <- sprintf("[data-for-btn=%s] .btn-loading-indicator", buttonId)
  doneEl <- sprintf("[data-for-btn=%s] .btn-done-indicator", buttonId)
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  shinyjs::disable(buttonId)
  shinyjs::show(selector = loadingEl)
  shinyjs::hide(selector = doneEl)
  shinyjs::hide(selector = errEl)
  on.exit({
    shinyjs::enable(buttonId)
    shinyjs::hide(selector = loadingEl)
  })
  
  # Try to run the code when the button is clicked and show an error message if
  # an error occurs or a success message if it completes
  tryCatch({
    value <- expr
    shinyjs::show(selector = doneEl, time = 30)
    shinyjs::delay(2000, shinyjs::hide(selector = doneEl, anim = TRUE, animType = "fade",
                                       time = 0.5))
    value
  }, error = function(err) { errorFunc(err, buttonId) })
}

errorFunc <- function(err, buttonId) {
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  errElMsg <- sprintf("[data-for-btn=%s] .btn-err-msg", buttonId)
  errMessage <- gsub("^ddpcr: (.*)", "\\1", err$message)
  shinyjs::html(html = errMessage, selector = errElMsg)
  shinyjs::show(selector = errEl, anim = TRUE, animType = "fade")
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
  if(!any(key == c("TSNE", "UMAP", "PCA", "DIFF")))
    stop("Reduction unknown")
  if(key == "PCA")
    return("PC")
  if(key == "DIFF")
    return("DC")
  if(key == "TSNE")
    return("tSNE")
  if(key == "UMAP")
    return("UMAP")
}

###### percentAbove #######
percentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}

######### dotPlot ############
dotPlot <- function(
  data,
  assay = NULL,
  features,
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
  group.by <- paste0(group.by,'_meta_data')
  if(!any(names(data$col.attrs) == group.by)){
    return(NULL)
  }
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  if(!any(names(data$col.attrs) == group.by) | !any(features%in%data$row.attrs$features[])){
    return(NULL)
  }
  index <- unlist(lapply(features, function(x){
    grep(x, data$row.attrs$features[], fixed=T)[1]
    }))
  data.features <- as.data.frame(data[['matrix']][,index])
  if(nrow(data.features)==0|ncol(data.features)==0) return(NULL)
  rownames(data.features) <- data$col.attrs$CellID[]
  colnames(data.features) <- features
  id <- data[[paste0("col_attrs/",group.by)]][drop=TRUE]
  if(length(unique(id))>50){
    showNotification(paste0('Notice: The dot plot will not load because the chosen grouping contains over 50 (',length(unique(id)),') unique groups.'), type = 'message')
    return(NULL)
  }
  data.features$id <- id
  
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- data[[paste0("col_attrs/",split.by)]][drop=TRUE]
    if (length(x = unique(x = splits)) > length(x = cols)) {
      stop("Not enought colors for the number of groups")
    }
    cols <- cols[1:length(x = unique(x = splits))]
    names(x = cols) <- unique(x = splits)
    data.features$id <- paste(data.features$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
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
  if (!is.null(x = split.by)) {
    splits.use <- vapply(
      X = strsplit(x = as.character(x = data.plot$id), split = '_'),
      FUN = '[[',
      FUN.VALUE = character(length = 1L),
      2
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        return(colorRampPalette(colors = c('grey', color))(20)[value])
      },
      color = cols[splits.use],
      value = avg.exp.scaled
    )
  }
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
  
  plot <- ggplot(data = data.plot, mapping = aes(x = data.plot$features.plot, y = data.plot$id)) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
    #scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(
      x = 'Features',
      y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
    )
  if (!is.null(x = split.by)) {
    plot <- plot + scale_color_identity()
  } else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (is.null(x = split.by)) {
    plot <- plot + guides(color = guide_colorbar(title = 'Average Expression'))
  }
  return(plot)
}

######## dimPlotlyOutput #######
dimPlotlyOutput <- function(assay.in, reduc.in, group.by, annot_panel = NULL, tmp_annotations = NULL, low.res, data){
  
  col.attrs <- names(data[[assay.in]]$col.attrs)

  group.by <- paste0(group.by,'_meta_data')
  # if(!grepl(paste0('_',tolower(assay.in),'_'),reduc.in)){
  #   return(NULL)
  # }
  if(!any(col.attrs == group.by)){
    return(NULL)
  }
  reduc <- col.attrs[grepl(reduc.in,col.attrs)]
  n <- length(reduc)
  
  if(n<2 || n>3){
    showNotification(ui = paste0('Error: Invalid number of dimensions for ', reduc.in, ': ',n),type = 'error')
  }

  plot.data <- data[[assay.in]]$get.attribute.df(attributes=c(reduc,group.by,'percent.mito_meta_data'))
  
  ann <- 50
  if(length(unique(plot.data[,ncol(plot.data)-1]))>ann){
    if (is.numeric(plot.data[,ncol(plot.data)-1])){
      showNotification(ui = paste0('Warning: The annotations you chose contain over ', ann, ' unique numeric annotations and may not be suitable for vizualization'),type = 'warning')
    }else{
      # Maybe check for this in the shiny server to not give the user the option to even select this... 
      showNotification(ui = paste0('Error: The annotations you chose contain over ', ann, ' unique non-numeric annotations and are not suitable for vizualization'),type = 'error')
      return(NULL)
    }
  }
  label_key <- reduc_key(key = toupper(sub("_.*","",reduc.in)))
  ax.x <- list(
    title = paste0(label_key,"_1"),
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'
  )
  ax.y <- list(
    title = paste0(label_key,"_2"),
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'
  )
  
  ax.z <- list(
    title = paste0(label_key,"_3"),
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'
  )
  #str(plot.data)
  #print(head(plot.data))
  if(!any(is.na(as.numeric(plot.data[,ncol(plot.data)-1]))) & length(unique(plot.data[,ncol(plot.data)-1,drop = TRUE]))<50){
    plot.data[,ncol(plot.data)-1] <- as.factor(as.numeric(plot.data[,ncol(plot.data)-1]))
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
                   sub('_meta_data','',group.by),": ", plot.data[,3],"\n</br>",
                   label_key,"_1: ", format(plot.data[,1],digits=3),"\n",
                   "</br>",label_key,"_2: ", format(plot.data[,2],digits=3), "\n",
                   "</br>percent.mt: ", format(plot.data[,4],digits=3), "%"), 
                 hovertemplate = paste0('<b>%{text}</b><extra></extra>')
    ) %>% layout(title = ifelse(test = annot_panel == 'cell_annotation_custom', yes = paste0(assay.in, " data coloured by custom annotations"), no = paste0(assay.in, " data coloured by ", sub('_meta_data','',group.by))) ,xaxis = ax.x, yaxis = ax.y, dragmode='lasso')
  }else{
    p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2], z = plot.data[,3],
                 type = 'scatter3d', mode = 'markers', key = ~rownames(plot.data), stroke = I('dimgrey'),
                 color = col, text =  ~paste0(
                   sub('_meta_data','',group.by),": ", plot.data[,4],"\n</br>",
                   label_key,"_1: ", format(plot.data[,1],digits=3),"\n",
                   "</br>",label_key,"_2: ", format(plot.data[,2],digits=3), "\n",
                   "</br>",label_key,"_3: ", format(plot.data[,3],digits=3), "\n",
                   "</br>percent.mt: ", format(plot.data[,5],digits=3), "%"), 
                 hovertemplate = paste0('<b>%{text}</b><extra></extra>')
    ) %>% layout(title = ifelse(test = annot_panel == 'cell_annotation_custom', yes = paste0(assay.in, " data coloured by custom annotations"), no = paste0(assay.in, " data coloured by ", sub('_meta_data','',group.by))) ,scene = list(xaxis = ax.x, yaxis = ax.y, zaxis = ax.z, dragmode='lasso'),legend = list(x = 100, y = 0.5))
  }
  if(low.res == 'yes'){
    return(p %>% toWebGL())
  }else{
    return(p)
  }
}

####### featurePlotlyOutput ##########
featurePlotlyOutput <- function(assay.in, reduc.in, group.by, feature.in, low.res, data){
  
  group.by <- paste0(group.by,'_meta_data')
  # if(!grepl(paste0('_',tolower(assay.in),'_'),reduc.in)){
  #   return(NULL)
  # }
  if(!any(names(data[[assay.in]]$col.attrs) == group.by)){
    return(NULL)
  }
  
  data.features <- as.data.frame(data[[assay.in]][['matrix']][,which(data[[assay.in]]$row.attrs$features[]%in%feature.in)])
  if(nrow(data.features)==0|ncol(data.features)==0) return(NULL)
  
  data.annot <- data[[assay.in]]$get.attribute.df(attributes=group.by)
  
  rownames(data.features) <- data[[assay.in]]$col.attrs$CellID[]
  colnames(data.features) <- feature.in
  
  dims <- names(data[[assay.in]]$col.attrs)[grepl(reduc.in,names(data[[assay.in]]$col.attrs))]
  
  plot.data <- data[[assay.in]]$get.attribute.df(attributes=dims)
  plot.data <- cbind(plot.data,data.annot)
  plot.data <- cbind(plot.data,data.features)
  
  label_key <- reduc_key(key = toupper(sub("_.*","",reduc.in)))
  
  ax.x <- list(
    title = paste0(label_key,"_1"),
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'
  )
  ax.y <- list(
    title = paste0(label_key,"_2"),
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'
  )
  ax.z <- list(
    title = paste0(label_key,"_3"),
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'
  )
  
  if(length(dims) == 2){
    p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2],
                 type = 'scatter', mode = 'markers', key = ~rownames(plot.data),
                 color = plot.data[,4],text =  ~paste0(
                   sub("_meta_data","",group.by),": ", plot.data[,3],"\n",
                   "</br>",'Expression',": ", format(plot.data[,4],digits=3),"\n",
                   "</br>",label_key,"_1: ", format(plot.data[,1],digits=3),"\n",
                   "</br>",label_key,"_2: ", format(plot.data[,2],digits=3), "\n"), 
                 hovertemplate = paste0('<b>%{text}</b>',
                                        '<extra></extra>')
    ) %>% layout(title = feature.in ,xaxis = ax.x, yaxis = ax.y, dragmode='lasso')
  }else{
    p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2], z = plot.data[,3],
                 type = 'scatter3d', mode = 'markers', key = ~rownames(plot.data),
                 color = plot.data[,5],text =  ~paste0(
                   sub("_meta_data","",group.by),": ", plot.data[,4],"\n",
                   "</br>",'Expression',": ", format(plot.data[,5],digits=3),"\n",
                   "</br>",label_key,"_1: ", format(plot.data[,1],digits=3),"\n",
                   "</br>",label_key,"_2: ", format(plot.data[,2],digits=3), "\n",
                   "</br>",label_key,"_3: ", format(plot.data[,3],digits=3), "\n"), 
                 hovertemplate = paste0('<b>%{text}</b>',
                                        '<extra></extra>')
    ) %>% layout(title = feature.in , scene = list(xaxis = ax.x, yaxis = ax.y, zaxis = ax.z), dragmode='lasso' ,legend = list(x = 100, y = 0.5))
  }
  if(low.res == 'yes'){
    return(p%>% toWebGL())
  }else{
    return(p)
  }
}

split_dot_plot <- function(data, 
                           features, 
                           group.by, 
                           split.by,
                           col.max = 2.5,
                           col.min = -2.5,
                           cols = c("blue","red"),
                           dot.min = 0,
                           dot.scale = 6,
                           scale.by = "radius",
                           scale.max = NA,
                           scale.min = NA){
  
  group.by <- paste0(group.by,'_meta_data')
  split.by <- paste0(split.by,'_meta_data')
  index <- unlist(lapply(features, function(x){
    grep(x, data$row.attrs$features[], fixed=T)[1]
    }))
  data.features <- as.data.frame(data[['matrix']][,index])
  if(nrow(data.features)==0|ncol(data.features)==0) return(NULL)
  rownames(data.features) <- data$col.attrs$CellID[]
  colnames(data.features)[1:length(features)] <- features
  
  id <- data[[paste0("col_attrs/",group.by)]][drop=TRUE]
  if(length(unique(id))>50){
    showNotification(paste0('Notice: The dot plot will not load because the chosen grouping contains over 50 (',length(unique(id)),') unique groups.'), type = 'message')
    return(NULL)
  }
  data.features$id <- id
  
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- data[[paste0("col_attrs/",split.by)]][drop=TRUE]
    #if (length(x = unique(x = splits)) > length(x = cols)) {
    #  stop("Not enought colors for the number of groups")
    #}
    cols <- cols[1:length(x = unique(x = splits))]
    names(x = cols) <- unique(x = splits)
    data.features$id <- paste(data.features$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
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
  
  data.plot$split <- sub(".*_","",data.plot$id)
  data.plot$id <- sub("_.*","",data.plot$id)
  
  data.plot <- data.plot
  data.plot$id <- as.factor(data.plot$id)
  data.plot$split <- as.factor(data.plot$split)
  
  data.plot$id <- reorder_levels(data.plot$id)
  #print('data.plot')
  #print(str(data.plot))
  
  labels <- data.frame(features.plot = rep(levels(data.plot$features.plot),length(levels(data.plot$split))))
  labels$id <- rep(levels(data.plot$id)[length(levels(data.plot$id))], nrow(labels))
  labels <- labels[order(labels$features.plot),]
  labels$split <- rep(levels(data.plot$split),nrow(labels)/length(levels(data.plot$split)))
  
  labels$features.plot <- as.factor(labels$features.plot)
  labels$id <- as.factor(labels$id)
  labels$split <- as.factor(labels$split)
  labels$y <- as.factor('ZZZZZZ') #length(levels(data.plot$id))+1
  #print('labels:')
  #print(labels)
  
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
  cat('\n\n')
  return(p)
}

