seurat_reduc_key_name <- function(object, reduction_name = NULL, key_as_name = TRUE, key_name = NULL){
  if(reduction_name==FALSE) stop('ERROR: you must specify the reduction name')
  if(key_as_name == FALSE && is.null(key_name)){
    stop('ERROR: If you do not want to set the key name to the reduction name, a key name must be specified')
  }else if(key_as_name == TRUE && !is.null(key_name)){
    key_name <- reduction_name
  }
  embeddings <- object@reductions[[reduction_name]]@cell.embeddings
  n <- ncol(embeddings)
  if(n<2 || n>3) stop(paste0('invalid number of reduced dimensions for visualization: ',n))
  colnames(embeddings) <- paste0(key_name,"_",1:n)
  object@reductions[[reduction_name]]@cell.embeddings <- embeddings
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
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  group.by <- paste0(group.by,'.viz')
  if(!any(names(data$col.attrs) == group.by) | !any(features%in%data$row.attrs$features[])){
    return(NULL)
  }
  data.features <- as.data.frame(data[['matrix']][,which(data$row.attrs$features[]%in%features)])
  if(nrow(data.features)==0|ncol(data.features)==0) return(NULL)
  rownames(data.features) <- data$col.attrs$CellID[]
  colnames(data.features)[1:length(features)] <- features
  data.features$id <- data[[paste0("col_attrs/",group.by)]][drop=TRUE]
  
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
  plot <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'id')) +
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
dimPlotlyOutput <- function(assay.in, reduc.in, group.by, data){
  
  col.attrs <- names(data[[assay.in]]$col.attrs)
  
  #print(paste(assay.in, reduc.in, group.by))
  group.by <- paste0(group.by,'_meta_data')
  if(!grepl(paste0('_',tolower(assay.in),'_'),reduc.in)){
    return(NULL)
  }
  if(!any(col.attrs == group.by)){
    return(NULL)
  }
  
  reduc <- col.attrs[grepl(reduc.in,col.attrs)]
  n <- length(reduc)
  if(n >2)

  if(group.by == "cluster"){
    plot.data <- data[[assay.in]]$get.attribute.df(attributes=c(reduc,tolower(paste0(assay.in,"_","clusters")),'percent.mt'))
  }
  else{
    plot.data <- data[[assay.in]]$get.attribute.df(attributes=c(reduc,group.by,'percent.mt'))
  }
  
  ann <- 50
  if(length(unique(plot.data[,3]))>ann){
    showNotification(ui = paste0('Warning: The annotations you chose contains over ', ann, 'unique annotations and may not be suitable for vizualization'),type = 'warning')
  }
  #print(head(plot.data))
  key <- reduc_key(key = toupper(sub("_.*","",reduc.in)))
  
  ax.x <- list(
    title = paste0(key,"_1"),
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'
  )
  ax.y <- list(
    title = paste0(key,"_2"),
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'
  )
  
  ax.z <- list(
    title = paste0(key,"_3"),
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'
  )
  
  #print(head(plot.data))
  if(n == 2){
    p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2], 
                 type = 'scatter', mode = 'markers', 
                 color = plot.data[,3], text =  ~paste0(
                   key,"_1: ", format(plot.data[,1],digits=3),"\n",
                   "</br>",key,"_2: ", format(plot.data[,2],digits=3), "\n",
                   "</br>percent.mt: ", format(plot.data[,4],digits=3), "%"), 
                 hovertemplate = paste0('<b>%{text}</b>')
    ) %>% layout(title = paste0(assay.in, " data coloured by ", sub('.viz','',group.by)) ,xaxis = ax.x, yaxis = ax.y)
  }else{
    p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2], z = plot.data[,3],
                 type = 'scatter3d', mode = 'markers', 
                 color = plot.data[,4], text =  ~paste0(
                   key,"_1: ", format(plot.data[,1],digits=3),"\n",
                   "</br>",key,"_2: ", format(plot.data[,2],digits=3), "\n",
                   "</br>",key,"_3: ", format(plot.data[,3],digits=3), "\n",
                   "</br>percent.mt: ", format(plot.data[,5],digits=3), "%"), 
                 hovertemplate = paste0('<b>%{text}</b>')
    ) %>% layout(title = paste0(assay.in, " data coloured by ", sub('.viz','',group.by)) ,scene = list(xaxis = ax.x, yaxis = ax.y, zaxis = ax.z),legend = list(x = 100, y = 0.5))
  }
  return(p)
}

####### featurePlotlyOutput ##########
featurePlotlyOutput <- function(assay.in, reduc.in, group.by, feature.in, data){
  
  group.by <- paste0(group.by,'.viz')
  if(!grepl(paste0('_',tolower(assay.in),'_'),reduc.in)){
    #print('1')
    return(NULL)
  }
  if(!any(names(data[[assay.in]]$col.attrs) == group.by)){
    print('3')
    return(NULL)
  }
  #print('4')
  
  
  data.features <- as.data.frame(data[[assay.in]][['matrix']][,which(data[[assay.in]]$row.attrs$features[]%in%feature.in)])
  if(nrow(data.features)==0|ncol(data.features)==0) return(NULL)
  if(group.by == 'cluster'){
    data.annot <- data[[assay.in]]$get.attribute.df(attributes=tolower(paste0(assay.in,"_clusters")))
  }else{
    data.annot <- data[[assay.in]]$get.attribute.df(attributes=group.by)
  }
  rownames(data.features) <- data[[assay.in]]$col.attrs$CellID[]
  colnames(data.features) <- feature.in
  
  n <- if(grepl("2d", reduc.in)){
    2
  }else{
    3
  }
  
  dims <- names(data[[assay.in]]$col.attrs)[grepl(reduc.in,names(data[[assay.in]]$col.attrs))]
  
  plot.data <- data[[assay.in]]$get.attribute.df(attributes=dims)
  plot.data <- cbind(plot.data,data.annot)
  plot.data <- cbind(plot.data,data.features)
  
  key <- reduc_key(key = toupper(sub("_.*","",reduc.in)))
  
  ax.x <- list(
    title = paste0(key,"_1"),
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'
  )
  ax.y <- list(
    title = paste0(key,"_2"),
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'
  )
  ax.z <- list(
    title = paste0(key,"_3"),
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = FALSE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='none'
  )
  
  if(n == 2){
    p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2],
                 type = 'scatter', mode = 'markers', 
                 color = plot.data[,4],text =  ~paste0(
                   key,"_1: ", format(plot.data[,1],digits=3),"\n",
                   "</br>",key,"_2: ", format(plot.data[,2],digits=3), "\n",
                   "</br>",group.by,": ", plot.data[,3]), 
                 hovertemplate = paste0('<b>%{text}</b>',
                                        '<extra></extra>')
    ) %>% layout(title = feature.in ,xaxis = ax.x, yaxis = ax.y)
  }else{
    p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2], z = plot.data[,3],
                 type = 'scatter3d', mode = 'markers', 
                 color = plot.data[,5],text =  ~paste0(
                   key,"_1: ", format(plot.data[,1],digits=3),"\n",
                   "</br>",key,"_2: ", format(plot.data[,2],digits=3), "\n",
                   "</br>",key,"_3: ", format(plot.data[,3],digits=3), "\n",
                   "</br>",group.by,": ", plot.data[,4]), 
                 hovertemplate = paste0('<b>%{text}</b>',
                                        '<extra></extra>')
    ) %>% layout(title = feature.in , scene = list(xaxis = ax.x, yaxis = ax.y, zaxis = ax.z), legend = list(x = 100, y = 0.5))
  }
  return(p)
}
