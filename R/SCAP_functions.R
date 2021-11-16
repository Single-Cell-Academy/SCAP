library("shiny")
library("cowplot")
library("ggplot2")
library("ggridges")
library("ggthemes")
library("gtools")
library("Matrix")
library("plotly")
library("Seurat")
library("viridis")
library("Nebulosa")
library("Seurat")
library("SeuratDisk")
library("reticulate")

COLORPAL_CONTINUOUS <- colorRampPalette(c("lightgrey", viridis(10)))

COLORPAL_DISCRETE <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941",
"#006FA6", "#A30059", "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", 
"#8FB0FF", "#997D87", "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", 
"#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", 
"#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", 
"#013349", "#00846F", "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", 
"#C2FF99", "#001E09", "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", 
"#7A87A1", "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", 
"#0086ED", "#886F4C", "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", 
"#FF913F", "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", 
"#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", 
"#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", 
"#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", 
"#324E72", "#6A3A4C", "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", 
"#806C66", "#222800", "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", 
"#1E0200", "#5B4E51", "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", 
"#7ED379", "#012C58", "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", 
"#837393", "#943A4D", "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", 
"#0AA3F7", "#E98176", "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", 
"#8D8546", "#9695C5", "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", 
"#55813B", "#E704C4", "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", 
"#643127", "#513A01", "#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", 
"#4C6001", "#9C6966", "#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", 
"#006C31", "#DDB6D0", "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", 
"#C6DC99", "#203B3C", "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", 
"#8BB400", "#797868", "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", 
"#CCB87C", "#B88183", "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", 
"#71BB8C", "#00B433", "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", 
"#BCB1E5", "#76912F", "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", 
"#A76F42", "#89412E", "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", 
"#EA8B66", "#958A9F", "#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", 
"#47675D", "#3A3F00", "#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", 
"#3C3E6E", "#D83D66", "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", 
"#866097", "#365D25", "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B")

use_raw <- function(h5ad){
  if(is.null(h5ad$raw)){
    return(FALSE) # no raw slot, don't use raw
  }else if(sum(h5ad$raw$X[1,]) %% 1 == 0){
    return(FALSE) # raw slot contains int counts, don't use raw
  }else{
    return(TRUE) # raw exists and doesn't contain int counts, use raw
  }
}

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

anndata_write_fix <- function(file){
  ad <- import('anndata')
  a = ad$read(file)
  py$tmp = a
  py_run_string("tmp.__dict__['_raw'].__dict__['_var'] = tmp.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})")
  if(!all(rownames(a$var) %in% rownames(a$raw$var))){
    py_run_string("tmp.__dict__['_raw'].__dict__['_var'] = tmp.__dict__['_raw'].__dict__['_var'].set_index('features', drop = False)")
    py_run_string("tmp.__dict__['_raw'].__dict__['_var'] = tmp.__dict__['_raw'].__dict__['_var'].rename_axis(index = None)")
  }
  a = py$tmp
  a$write(file)
}

loom_to_h5ad <- function(file_1, file_2){
  message("Converting to h5ad")
  ad <- import('anndata')
  a = ad$read_loom(file_1)
  if(length(a$obsm_keys()) == 0){
    reductions <- a$obs_keys()[grep("_reduction$", a$obs_keys())]
    if(length(reductions) == 0){
      print('No dimensional reductions found in data')
    }else{
      reduction_names <- unique(sub("_\\d_reduction$", "", reductions))
      reduction_lst <- lapply(reduction_names, function(x) reductions[grep(x, reductions)])
      names(reduction_lst) <- reduction_names
      a = ad$read_loom(file_1, obsm_names = reduction_lst)
    }
  }
  if(!grepl("\\.h5ad", file_2)){
    file_2 <- paste0(file_2, ".h5ad")
  }
  a$write(file_2)
}

rds_to_h5ad <- function(file_1, file_2){
  obj <- readRDS(file_1)
  obj <- try_seurat_update(obj)
  #convert factors to characters
  idx <- which(unlist(lapply(obj@meta.data, class)) == 'factor')
  if(length(idx)>0){
    for(i in idx){
      obj@meta.data[,i] <- as.character(obj@meta.data[,i,drop=TRUE])
    }
  }
  assays <- names(obj@assays)
  message("Converting to h5seurat")
  for(assay in assays){
    print(length(obj[[assay]]@var.features))
    if(length(obj[[assay]]@var.features) == 0){
      obj <- FindVariableFeatures(obj)
    }
  }
  if(!grepl("\\.h5ad", file_2)){
    file_2 <- paste0(file_2, ".h5ad")
  }
  SaveH5Seurat(obj, filename = gsub("\\.h5ad$", ".h5Seurat", file_2), overwrite = TRUE)
  message("Converting to h5ad")
  for(assay in assays){
    message(paste0("Converting ", assay, "..."))
    Convert(gsub("\\.h5ad$", ".h5Seurat", file_2), dest = "h5ad", assay = assay, overwrite = TRUE)
    system(paste0("mv ", file_2, " ", paste0(gsub("\\.h5ad", "", file_2), "_", assay, '.h5ad')))
    anndata_write_fix(paste0(gsub("\\.h5ad", "", file_2), "_", assay, '.h5ad'))
  }
}

h5ad_to_rds <- function(file_1, file_2){
  message("Converting to h5seurat")
  Convert(file_1, dest = "h5seurat", overwrite = TRUE)
  message("Converting to rds")
  seur <- LoadH5Seurat(sub("\\.h5ad$", '.h5seurat', file_1))
  if(!grepl("\\.rds", file_2)){
    file_2 <- paste0(file_2, ".rds")
  }
  saveRDS(seur, file_2)
}

SCAP_Convert <- function(from, to, file_1, file_2){
  if(from == 'loom' & to == 'h5ad'){
    loom_to_h5ad(file_1, file_2)
  }else if(from == 'rds' & to == 'h5ad'){
    rds_to_h5ad(file_1, file_2)
  }else if(from == 'h5ad' & to == 'rds'){
    h5ad_to_rds(file_1, file_2)
  }else{
    message('error')
  }
}

# #------ required to fix anndata.write ValueError: '_index' is a reserved name for dataframe columns. Above error raised while writing key 'raw/var' of <class 'h5py._hl.files.File'> from /. ------#
# cat(file = stderr(), "raw/var correction...\n")
# a = anndata$read(out_path)
# py$tmp = a
# py_run_string("tmp.__dict__['_raw'].__dict__['_var'] = tmp.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})")
# a = py$tmp
# a$write(out_path)
# #--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

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
dimPlotlyOutput <- function(assay.in, reduc.in, group.by, annot_panel = NULL, tmp_annotations = NULL, low.res, hide.legend){

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
                 color = col, colors = COLORPAL_DISCRETE[1:length(unique(col))], text =  ~paste0(
                   names(group.by)[1],": ", plot.data[,3],"\n</br>",
                   colnames(reduc.in)[1],": ", format(plot.data[,1],digits=3),"\n",
                   "</br>",colnames(reduc.in)[2],": ", format(plot.data[,2],digits=3)),
                 hovertemplate = paste0('<b>%{text}</b><extra></extra>')
    ) %>% layout(title = ifelse(test = annot_panel == 'cell_annotation_custom', yes = paste0(assay.in, " data coloured by custom annotations"), no = paste0(assay.in, " data coloured by ", names(group.by))) ,xaxis = ax.x, yaxis = ax.y, dragmode='lasso')
  }else{
    p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2], z = plot.data[,3],
                 type = 'scatter3d', mode = 'markers', key = ~rownames(plot.data), stroke = I('dimgrey'),
                 color = col, colors = COLORPAL_DISCRETE[1:length(unique(col))], text =  ~paste0(
                   names(group.by)[1],": ", plot.data[,4],"\n</br>",
                   colnames(reduc.in)[1],": ", format(plot.data[,1],digits=3),"\n",
                   "</br>",colnames(reduc.in)[2],": ", format(plot.data[,2],digits=3), "\n",
                   "</br>",colnames(reduc.in)[3],": ", format(plot.data[,3],digits=3)), 
                 hovertemplate = paste0('<b>%{text}</b><extra></extra>')
    ) %>% layout(title = ifelse(test = annot_panel == 'cell_annotation_custom', yes = paste0(assay.in, " data coloured by custom annotations"), no = paste0(assay.in, " data coloured by ", names(group.by))) ,scene = list(xaxis = ax.x, yaxis = ax.y, zaxis = ax.z, dragmode='lasso'),legend = list(x = 100, y = 0.5))
  }

  if(low.res == 'yes'){
    p <- p %>% toWebGL()
  }
  
  if(hide.legend == 'Yes'){
    p <- p %>% layout(showlegend = FALSE)
  }
  
  return(p)
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
                 color = plot.data[,4], colors = COLORPAL_CONTINUOUS(100), opacity = 0.6, text =  ~paste0(
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
                 color = plot.data[,5], colors = COLORPAL_CONTINUOUS(100), opacity = 0.6, text =  ~paste0(
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
    return(plot_density(data_seurat, features = colnames(feature.in), size = 1) + theme_cowplot())
  }else{
    p <- plot_density(data_seurat, features = colnames(feature.in), joint = TRUE, combine = FALSE, size = 1) 
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
  labels$avg.exp <- 0
  labels$pct.exp <- NA
  labels$avg.exp.scaled <- data.plot$avg.exp.scaled[1]
  labels$id <- as.factor('ZZZZZZ')
  labels <- labels[order(labels$features.plot),]
  labels$split <- rep(levels(data.plot$split),nrow(labels)/length(levels(data.plot$split)))
  
  labels$features.plot <- as.factor(labels$features.plot)
  labels$id <- as.factor(labels$id)
  labels$split <- as.factor(labels$split)
  
  data.plot <- rbind(data.plot, labels[,c(2,3,1,4:ncol(labels))])
  tmp <<- data.plot
  tmp_labs <<- labels
  p <- ggplot(data = data.plot, mapping = aes(x = features.plot, y = id)) + 
    geom_point(aes(color = avg.exp.scaled, size = pct.exp, group = split), position = position_dodge(1))  +
    geom_text(data = data.plot[which(data.plot$id == "ZZZZZZ"),], aes(label = split, group = split), position = position_dodge(width = 1), size=5) + 
    scale_y_discrete(breaks = levels(data.plot$id)[1:(length(levels(data.plot$id))-1)]) +
    scale_color_gradient(low = 'blue', high = 'red') + 
    xlab('Features') +
    ylab('Identity') +
    theme_base()
  tmp_p <<- p
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
  tmp_p_2 <<- p
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

  #cat(file = stderr(), "data$expression: ", summary(data$expression), "\n")

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
  geom_density_ridges(alpha=0.6, bandwidth=0.5,
                      quantile_lines = TRUE, quantiles = 2) +
  xlab("Expression") +
  ylab("Identity") + 
  facet_grid(. ~ feature, scales = (if (same.y.lims) 'fixed' else 'free')) +
  FacetTheme(
    panel.spacing = unit(0, 'lines'),
    panel.background = element_rect(fill = NA, color = "black"),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y.right = element_text(angle = 0)
  ) + NoLegend()

  return(plot)
}

