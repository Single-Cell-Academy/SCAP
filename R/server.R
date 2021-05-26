library("cowplot")
library("dplyr")
library("ggplot2")
library("ggthemes")
library("gtools")
library("Matrix")
library("MODIS")
library("plotly")
library("rjson")
library("shiny")
library("shinyFiles")
library("data.table")
library("scibet")
library("readr")
library("reactable")
library("reticulate")
library("shinyjs")
library("presto")
library("bbplot")

reticulate::use_virtualenv("../renv/python/virtualenvs/renv-python-3.8.5/")

#### Variables that persist across sessions
## Read in table with datasets available for SciBet
datasets_scibet <- fread("../meta/SciBet_reference_list.tsv")
## Source functions
source("SCAP_functions.R")
source_python("../Python/rank_genes_groups_df.py")

anndata <- import('anndata')
scanpy <- import('scanpy')

init <- 0 # flag for autosave

server <- function(input, output, session){
  session$onSessionEnded(stopApp)
  
  options(shiny.maxRequestSize=500*1024^2)

  rvalues <- reactiveValues(tmp_annotations = NULL, cells = NULL, order = NULL, features = NULL, obs = NULL, obs_cat = NULL, reductions = NULL, cell_ids = NULL, h5ad = NULL, path_to_data = NULL,
                            raw_dtype = NULL)
  rvalues_mod <- reactiveValues(tmp_annotations = NULL, cells = NULL, order = NULL, features = NULL, obs = NULL, obs_cat = NULL, reductions = NULL, cell_ids = NULL, h5ad = NULL, path_to_data = NULL,
                                raw_dtype = NULL)

  ## Determine folders for ShinyDir button
  volumes <- c("FTP" = "/ftp", Home = fs::path_home())
  
  ## GenAP2 logo
  output$genap_logo <- renderImage({
    # Return a list containing the filename
    list(src = "./img/GenAP_powered_reg.png",
     contentType = 'image/png',
     width = "100%",
     height = "100%",
     alt = "This is alternate text")
    }, deleteFile = FALSE)

  ## File directory
  shinyFileChoose(input, "h5ad_in", roots = volumes, session = session)
  
  # connect chosen .h5ad file
  observeEvent(input$h5ad_in, {
    path <- parseFilePaths(selection = input$h5ad_in, roots = volumes)$datapath
    if(is.integer(path[1]) || identical(path, character(0)) || identical(path, character(0))) return(NULL)
    h5ad_files <- path#paste0(path,"/",list.files(path))
    assays <- sub(".h5ad","",sub(paste0(".*/"),"",h5ad_files))
    data <- list()
    ## Iterate over all assays and connect to h5ad objects
    for(i in 1:length(assays)){
      data[[i]] <- tryCatch({
        anndata$read(h5ad_files[i])
      },
      error = function(e){
        showModal(modalDialog(p(paste0("An error occured trying to connect to ", h5ad_files[i])), title = "Error connecting to h5ad file."), session = getDefaultReactiveDomain())
        return(NULL)
      })
    }
    if(is.null(data)) return(NULL)
    if(length(data) != length(assays)) return(NULL)
    if(length(unlist(lapply(data, function(x){x}))) != length(assays)) return(NULL)
    names(data) <- assays
    ## Check if RAW Anndata object is present or not. If not present, use the main object
    if(is.null(data[[1]]$raw)){ 
      rvalues$features <- rownames(data[[1]]$var)
    }else{
      test_gene_name <- rownames(data[[1]]$var)[1]
      if(test_gene_name %in% rownames(data[[1]]$raw$var)){ # check if rownames are numbers or gene names
        rvalues$features <- rownames(data[[1]]$raw$var)
      }else if("features" %in% colnames(data[[1]]$raw$var)){ ## Check if there is a column named features in raw
        rvalues$features <- data[[1]]$raw$var$features
      }else if(test_gene_name %in% data[[1]]$raw$var[,1]){ # otherwise, check if the first column contains rownames
        rvalues$features <- data[[1]]$raw$var[,1]
      }
    }
    rvalues$obs <- data[[1]]$obs_keys()
    ## Determine type of annotation and create a layer to annotate for easy usage later on
    rvalues$obs_cat <- check_if_obs_cat(obs_df = data[[1]]$obs) ## Function to check if an observation is categorical or numeric
    reductions <- data[[1]]$obsm$as_dict()
    reduction_keys <- data[[1]]$obsm_keys()
    r_names <- rownames(data[[1]]$obs)
    for(i in 1:length(reductions)){
      reductions[[i]] <- as.data.frame(reductions[[i]])
      colnames(reductions[[i]]) <- paste0(reduction_keys[i], "_", 1:ncol(reductions[[i]]))
      rownames(reductions[[i]]) <- r_names
    }
    names(reductions) <- reduction_keys
    rvalues$reductions <- reductions
    rvalues$cell_ids <- rownames(data[[1]]$obs)
    rvalues$h5ad <- data
    rvalues$path_to_data <- h5ad_files
    ## Determine what data is likely stored in .raw
    if(is.null(data[[1]]$raw)){ ## Check if raw exists
      rvalues$raw_dtype <- "NULL"
    }else if(sum(rvalues$h5ad[[1]]$raw$X[1,]) %% 1 == 0){ ## Check whether raw contains un-normalized data or normalized data
      rvalues$raw_dtype <- "counts"
    }else{ ## Only if the other two conditions fail, use raw values to calculate differential expression
      rvalues$raw_dtype <- "normalized"
    }
    
    init <<- 0
  })

  # observe({ # auto save h5ad file(s)
  #   req(rvalues$h5ad)
  #   invalidateLater(120000) # 2 min
  #   if(init>0){
  #     #tryCatch(
  #     #  { 
  #         cat(file = stderr(), paste0(rvalues$path_to_data, "\n"))
  #         showNotification("Saving...", duration = NULL, id = 'auto_save')
  #         for(i in 1:length(rvalues$path_to_data)){
  #           rvalues$h5ad[[i]]$write(filename = rvalues$path_to_data[i])
  #         }
  #         removeNotification(id = 'auto_save')
  #      # },
  #      # error = function(e)
  #      # {
  #         #cat(file = stderr(), unlist(e))
  #     #    showModal(modalDialog(p(paste0("An error occured trying to write to ", rvalues$path_to_data[i], ": ", unlist(e))), title = "Error writing to h5ad file."), session = getDefaultReactiveDomain())
  #     #  }
  #    # )
  #   }
  #   init <<- init + 1
  # })

  ###==============// MAIN TAB //==============####
  
  #-- select input for Assay --#
  output$assay_1 <- renderUI({
    req(rvalues$h5ad)
    selectInput('assay_1', "Select Assay", 
      choices = names(rvalues$h5ad), 
      selected = ifelse(any(names(rvalues$h5ad)=="RNA"),
        yes = "RNA",
        no = names(rvalues$h5ad)[1])
      )
  })

  #-- Display name of data --#
  output$data_used <- renderText({
    if(!is.null(input$assay_1) & !is.na(input$h5ad_in[[1]][2])){
      path <- sub(".*\\/", "", parseFilePaths(selection = input$h5ad_in, roots = volumes)$datapath)
      return(paste0("Chosen data: ", path))
    }else{
      return(paste0(""))
    }
  })
  
  #-- select how to group cells --#
  output$grouping_1 <- renderUI({
    req(input$assay_1, rvalues$obs)
    assay <- input$assay_1
    if(any(grepl("seurat_clusters", rvalues$obs, ignore.case = TRUE))){
      sel <- rvalues$obs[grep("seurat_clusters", rvalues$obs, , ignore.case = TRUE)]
    }else if(any(grepl(paste0(tolower(assay),"_clusters"), rvalues$obs, ignore.case = TRUE))){
      sel <- rvalues$obs[grep(paste0(tolower(assay),"_clusters"), rvalues$obs, , ignore.case = TRUE)]
    }else if(any(grepl("louvain", rvalues$obs, ignore.case = TRUE))){
      sel <- rvalues$obs[grep("louvain", rvalues$obs, ignore.case = TRUE)]
    }else{
      sel <- rvalues$obs[1]
    }
    selectInput(inputId = 'grouping_1', label = 'Group By', choices = rvalues$obs, selected = sel, multiple = FALSE)
  })
  
  #-- select the clustering method --#
  output$reduction_1 <- renderUI({
    req(input$assay_1)
    assay <- input$assay_1
    options <- names(rvalues$reductions)
    sel <- if(any(grepl("umap", options, ignore.case = TRUE))){
      options[grepl("umap", options, ignore.case = TRUE)][1]
    }else if(any(grepl("tsne", options, ignore.case = TRUE))){
      options[grepl("tsne", options, ignore.case = TRUE)][1]
    }else{
      options[1]
    }
    selectInput(inputId = 'reduction_1', 'Choose Clustering Method', choices = as.list(options), selected = sel)
  })
  
  #-- select the features for the feature plot --#
  output$featureplot_1_feature_select <- renderUI({
    req(input$assay_1)
    assay <- input$assay_1
    selectInput(inputId = 'featureplot_1_feature_select', 
                label = 'Select a Feature to Visualize on the Feature Plot', 
                choices = rvalues$features, 
                selected = rvalues$features[1], 
                multiple = ifelse(input$nebulosa_on == "yes", TRUE, FALSE))
  })
  
  #-- select the featutres for the dot plot --#
  output$dotplot_1_feature_select<- renderUI({
    req(input$assay_1)
    assay <- input$assay_1

    selectInput(inputId = 'dotplot_1_feature_select', 
                label = 'Choose Features to Visualize on Dot Plot', 
                choices = rvalues$features, 
                selected = rvalues$features[1:2], 
                multiple = TRUE)
  })
  
  #-- visualize dotplot as a split dotplot --#
  output$do_split <- renderUI({
    req(input$assay_1)
    radioButtons(inputId = 'do_split', label = 'Split dot plot', choices = c('yes', 'no'), selected = 'no', inline = TRUE)
  })
  
  #-- select how to split the dot plot --#
  output$split_by <- renderUI({
    req(input$assay_1, rvalues$obs)
    choices <- unlist(lapply(rvalues$obs[which(rvalues$obs_cat)], function(x){
      if(length(unique(rvalues$h5ad[[1]]$obs[x][,,drop=TRUE]))>1){
        return(x)
      }
    }))
    if(is.null(choices)){
      return("No conditions to split by.")
    }else{
      return(radioButtons(inputId = 'split_by', label = 'Choose how to split the data', choices = choices, inline = FALSE))
    }
  })
  
  ## output variable to hold type of user selected grouping for main panel
  output$grouping_1_type <- reactive({
    req(input$grouping_1)
    cat <- rvalues$obs_cat[which(rvalues$obs == input$grouping_1)]
    if(cat){"yes"}else{"no"}
  })
  outputOptions(output, "grouping_1_type", suspendWhenHidden = FALSE)
  
  #-- dimensional reduction plot coloured by cell groups --#
  output$dimplot_1 <- renderPlotly({
    req(input$grouping_1, input$reduction_1)
    group.by <- list(rvalues$h5ad[[1]]$obs[input$grouping_1][,,drop=TRUE])
    names(group.by) <- input$grouping_1
    names(group.by[[1]]) <- rvalues$cell_ids
    cat <- rvalues$obs_cat[which(rvalues$obs == input$grouping_1)]
    if(cat){
      dimPlotlyOutput(assay.in = input$assay_1, 
                    reduc.in = rvalues$reductions[[input$reduction_1]], 
                    group.by = group.by, 
                    annot_panel = "", 
                    low.res = 'yes')
    }else{
      feature.in <- rvalues$h5ad[[1]]$obs[input$grouping_1][,,drop=FALSE]
      colnames(feature.in) <- input$grouping_1
      rownames(feature.in) <- rownames(rvalues$reductions[[input$reduction_1]])
      featurePlotlyOutput(assay.in = input$assay_1,
                          reduc.in = rvalues$reductions[[input$reduction_1]],
                          group.by = group.by,
                          feature.in = feature.in,
                          low.res = 'yes')
    }
  })
  
  #-- dimensional reduction plot coloured by feature expression --#
  output$featureplot_1 <- renderPlotly({
    req(input$grouping_1, input$reduction_1, input$featureplot_1_feature_select)
    group.by <- list(rvalues$h5ad[[1]]$obs[input$grouping_1][,,drop=TRUE])
    names(group.by) <- input$grouping_1
    names(group.by[[1]]) <- rvalues$cell_ids
    if(rvalues$raw_dtype == "NULL" | rvalues$raw_dtype == "counts"){
      feature.in <- as.data.frame(rvalues$h5ad[[1]]$X[,match(input$featureplot_1_feature_select, rvalues$features)])
    }else if(rvalues$raw_dtype == "normalized"){
      feature.in <- as.data.frame(rvalues$h5ad[[1]]$raw$X[,match(input$featureplot_1_feature_select, rvalues$features)])
    }
    colnames(feature.in) <- input$featureplot_1_feature_select
    rownames(feature.in) <- rownames(rvalues$reductions[[input$reduction_1]])
    if(input$nebulosa_on == 'no'){
      featurePlotlyOutput(assay.in = input$assay_1,
                          reduc.in = rvalues$reductions[[input$reduction_1]],
                          group.by = group.by,
                          feature.in = feature.in,
                          low.res = 'yes')
    }else{
      featurePlotlyOutput_nebulosa(assay.in = input$assay_1,
                                reduc.in = rvalues$reductions[[input$reduction_1]],
                                group.by = as.data.frame(group.by),
                                feature.in = feature.in,
                                low.res = 'yes')
    }
  })
  
  #-- dot plot for selected feature expression --#
  dot_plot <- reactive({
    req(input$dotplot_1_feature_select)

    assay <- input$assay_1

    ax.x <- list(
    title = "Features",
    zeroline = FALSE,
    showline = TRUE,
    showticklabels = TRUE,
    showgrid = FALSE,
    mirror=TRUE,
    ticks='outside',
    tickangle = -45
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
    
    if(rvalues$raw_dtype == "NULL" | rvalues$raw_dtype == "counts"){
      data.features <- as.data.frame(rvalues$h5ad[[1]]$X[,match(input$dotplot_1_feature_select, rvalues$features)])
    }else if(rvalues$raw_dtype == "normalized"){
      data.features <- as.data.frame(rvalues$h5ad[[1]]$raw$X[,match(input$dotplot_1_feature_select, rvalues$features)])
    }
    colnames(data.features) <- input$dotplot_1_feature_select
    rownames(data.features) <- rownames(rvalues$reductions)
    data.features$id <- rvalues$h5ad[[1]]$obs[input$grouping_1][,,drop=TRUE]

    if(input$do_split == 'yes' & !is.null(input$split_by)){
      splits = list(rvalues$h5ad[[1]]$obs[input$split_by][,,drop=TRUE])
      names(splits) <- input$split_by
      p <- split_dot_plot(data.features = data.features, features = input$dotplot_1_feature_select, assay = assay, split.by = splits)
      if(is.null(p)) return(NULL)
      if(identical(class(p),"shiny.tag")) return(NULL)
      return(p)
    }else if(input$do_split == 'no'){
      p <- dotPlot(data.features = data.features, features = input$dotplot_1_feature_select, assay = assay)
      if(is.null(p)){
        return(NULL)
      }
      if(identical(class(p),"shiny.tag")) return(NULL)
      p <- p + guides(color = guide_colourbar(order = 1, title = "Average Expression"), size = guide_legend(order = 2, title = "Percent Expressed")) + theme_few()
      return(p)
    }else{
      return(-1)
    }
  })
  
  output$dotplot_1 <- renderUI({
    req(dot_plot())
    if(input$do_split == 'yes' && is.null(input$split_by)){
      output$choose_splits <- renderText({
        return('Please choose how to split the data')
      })
      return(h2(textOutput('choose_splits')))
    }else{
      output$plot <- renderPlot({
        return(dot_plot())
      })
      return(plotOutput('plot'))
    }
  })
  
  
  ###==============// ANNOTATION TAB //==============####
  
  #-- select input for Assay --#
  output$assay_2 <- renderUI({
    req(rvalues$h5ad)
    selectInput('assay_2', "Select Assay", 
      choices = names(rvalues$h5ad), 
      selected = ifelse(any(names(rvalues$h5ad)=="RNA"),
        yes = "RNA",
        no = names(rvalues$h5ad)[1]))
  })
  
  to_listen <- reactive({
    list(input$assay_2, rvalues$h5ad)
  })

  #-- initialize/reset tmp_annotations --#
  observeEvent(to_listen(),{
    req(input$assay_2, rvalues$h5ad)
    rvalues$tmp_annotations <- rep("unlabled", rvalues$h5ad[[1]]$n_obs)
    names(rvalues$tmp_annotations) <- rownames(rvalues$h5ad[[1]]$obs)
  })
  
  #-- select how to group cells --#
  output$grouping_2 <- renderUI({
    req(input$assay_2, rvalues$obs)
    assay <- input$assay_2
    options <- rvalues$obs[which(rvalues$obs_cat)] # only show categorical metadata
    if(any(grepl("seurat_clusters", options, ignore.case = FALSE))){
      sel <- rvalues$obs[grep("seurat_clusters", options)]
    }else if(any(grepl(paste0(tolower(assay),"_clusters"), options, ignore.case = TRUE))){
      sel <- paste0(tolower(assay),"_clusters")
    }else{
      sel <- options[1]
    }
    selectInput(inputId = 'grouping_2', label = 'Group By', choices = options, selected = sel, multiple = FALSE)
  })
  
  #-- select the clustering method --#
  output$reduction_2 <- renderUI({
    req(input$assay_2)
    assay <- input$assay_2
    options <- names(rvalues$reductions)
    sel <- if(any(grepl("umap", options, ignore.case = TRUE))){
      options[grepl("umap", options, ignore.case = TRUE)][1]
    }else if(any(grepl("tsne", options, ignore.case = TRUE))){
      options[grepl("tsne", options, ignore.case = TRUE)][1]
    }else{
      options[1]
    }
    selectInput(inputId = 'reduction_2', 'Choose Clustering Method', choices = as.list(options), selected = sel)
  })
  
  #-- select the features for the feature plot --#
  output$featureplot_2_feature_select <- renderUI({
    req(input$assay_2)
    assay <- input$assay_2
    selectInput(inputId = 'featureplot_2_feature_select', 
                label = 'Select a Feature to Visualize on the Feature Plot', 
                choices = rvalues$features, 
                selected = rvalues$features[1], 
                multiple = ifelse(input$nebulosa_on == "yes", TRUE, FALSE))
  })
  
  #-- display cells that were selected in a table --#
  output$selected_cells <- renderTable({
    req(input$assay_2)
    selected_cells <- as.data.frame(event_data("plotly_selected")$key, stringsAsFactors = FALSE)
    if (is.null(selected_cells) | ncol(selected_cells) == 0 | nrow(selected_cells) == 0){
      rvalues$cells <- NULL
      return("Click and drag on either plot to select cells for marker identification (i.e., select/lasso) (double-click on either plot to clear selection)")
    }else if(ncol(selected_cells)>1){
      selected_cells <- as.data.frame(as.character(selected_cells[1,]),stringsAsFactors = F)
    }
    colnames(selected_cells) <- ""
    rvalues$cells <- selected_cells
    selected_cells
  })
  
  #-- dimensional reduction plot coloured by cell groups --#
  output$dimplot_2 <- renderPlotly({
    req(input$grouping_2, input$reduction_2)
    group.by <- list(rvalues$h5ad[[1]]$obs[input$grouping_2][,,drop=TRUE])
    names(group.by) <- input$grouping_2
    names(group.by[[1]]) <- rvalues$cell_ids
    dimPlotlyOutput(assay.in = input$assay_2, 
                    reduc.in = rvalues$reductions[[input$reduction_2]], 
                    group.by = group.by, 
                    annot_panel = input$annot_panel, 
                    tmp_annotations = rvalues$tmp_annotations, 
                    low.res = 'yes')
  })
  
  #-- dimensional reduction plot coloured by feature expression --#
  output$featureplot_2 <- renderPlotly({
    req(input$grouping_2, input$reduction_2, input$featureplot_2_feature_select)
    group.by <- list(rvalues$h5ad[[1]]$obs[input$grouping_2][,,drop=TRUE])
    names(group.by) <- input$grouping_2
    names(group.by[[1]]) <- rvalues$cell_ids
    if(rvalues$raw_dtype == "NULL" | rvalues$raw_dtype == "counts"){
      feature.in <- as.data.frame(rvalues$h5ad[[1]]$X[,match(input$featureplot_2_feature_select, rvalues$features)])
    }else if(rvalues$raw_dtype == "normalized"){
      feature.in <- as.data.frame(rvalues$h5ad[[1]]$raw$X[,match(input$featureplot_2_feature_select, rvalues$features)])
    }
    colnames(feature.in) <- input$featureplot_2_feature_select
    rownames(feature.in) <- rownames(rvalues$reductions[[input$reduction_2]])
    featurePlotlyOutput(assay.in = input$assay_2, 
                        reduc.in = rvalues$reductions[[input$reduction_2]], 
                        group.by = group.by, 
                        feature.in = feature.in, 
                        low.res = 'yes')
  })
  
  #-- Find the markers for the 1) selected cells compared to all other cells or
  #-- 2) the markers that define each group. Depending on if cells are selected or not
  observeEvent(input$find.markers, ignoreNULL = FALSE, ignoreInit = FALSE, handlerExpr = {
    if(input$find.markers == 0){ # Don't run before button click
      output$markers <- NULL
    }else{
      output$markers <- DT::renderDataTable({
        req(input$assay_2, rvalues$obs)
        cells <- isolate(rvalues$cells) # cell ids selected from scatter plots. Isolated so it is only triggered on button click

        y <- as.character(rvalues$h5ad[[1]]$obs[input$grouping_2][,,drop = TRUE]) # character vector of cell identities from selected meta data slot

        if(!is.null(cells)){ # if cells were selected on the scatter plots
          cols <- match(cells[,1], rvalues$cell_ids) # get index of selected cell ids from the main object
          y[cols] <- 'Selected' # rename the identity of the selected cells
        }

        if(length(unique(y))==1){ # need more than 1 cell identity for DE analysis
          showModal(modalDialog(p("Must have at least 2 groups defined to use Find Markers."), title = "Warning"), session = getDefaultReactiveDomain())
          return(NULL)
        }

        rvalues$h5ad[[1]]$obs['scap_find_markers_groups'] <- reorder_levels(y) # add cell identities to main object so scanpy methods can be used

        scanpy$tl$rank_genes_groups(rvalues$h5ad[[1]], groupby = 'scap_find_markers_groups', use_raw = TRUE, method = 'wilcoxon') # find DEGs

        t <- rank_genes_groups_df(rvalues$h5ad[[1]]) # Create a dataframe of DE results
        t$group <- py_to_r(attr(t, which='pandas.index')$tolist()) # some crpytic code that's required so an error doesn't occur
        colnames(t) <- c('feature', 'score', 'pval', 'pval_adj', 'logFC', 'group')
        t <- t[,c(1, ncol(t), 2:(ncol(t)-1))]
        if(!is.null(cells)){ # if cells were selected on the scatter plots, only show these cells.
          t <- t[which(t$group == "Selected"),]
        }
        return(t %>% arrange(desc(logFC)) %>% DT::datatable(filter = 'top') %>% DT::formatRound(columns=c('score', 'pval', 'pval_adj', 'logFC'), digits=3))
        #return(t %>% arrange(desc(Specificity)) %>% DT::datatable(filter = 'top') %>% DT::formatRound(columns=c("avgExpr", "logFC", "pval", "padj", "pct_in", "pct_out", "Specificity"), digits=3))
      })
    }
  })

  #-- display selected grouping name --#
  output$annotation_title <- renderText({
    req(input$grouping_2, rvalues$obs)
    return(paste0(input$grouping_2))
  })
  
  #-- display text boxes for user to enter new IDs for the meta data --#
  output$annotations <- renderUI({
    req(input$assay_2, input$grouping_2)
    names <- mixedsort(unique(rvalues$h5ad[[1]]$obs[input$grouping_2][,,drop = TRUE]))
    id <- "my_annot"
    lapply(1:length(names), function(x){
      textInput(inputId = paste0(names[x],"_",id), label = names[x], placeholder = names[x])
    })
  })
  
  #-- add user defined annotations to selected cells --#
  observeEvent(input$add_to_tmp,{
    if(is.null(rvalues$cells)){
      showNotification("Please select cells from either plot first", type = 'warning')
    }else if(is.null(input$custom_name) | input$custom_name == ""){
      showNotification("Please enter an annotation name", type = 'warning')
    }else{
      rvalues$tmp_annotations[which(names(rvalues$tmp_annotations)%in%rvalues$cells[,1])] <- input$custom_name
    }
  })

  #-- store user defined annotations --#
  observeEvent(input$set_cell_types, ignoreNULL = FALSE, ignoreInit = FALSE, handlerExpr = {
    req(rvalues$obs)
    group.by <- paste0(input$grouping_2)
    id <- "my_annot"
    if(is.null(input$new_scheme_name) || input$new_scheme_name == "" || input$new_scheme_name == " "){
      showNotification("Please name the annotation scheme", type = 'warning')
    }else if(input$new_scheme_name %in% rvalues$obs){
      showModal(modalDialog(p(paste0("ERROR ", input$new_scheme_name, " is already in the metadata. Please choose a unique name.")), title = "Error adding metadata"), session = getDefaultReactiveDomain())
    }else{
      if(input$annot_panel == 'cell_annotation_cluster'){
        names <- rvalues$h5ad[[1]]$obs[group.by][,,drop = TRUE]
        new <- names
        names <- unique(names)
        for(i in 1:length(names)){
          new[which(new == names[i])] <- input[[paste0(names[i],"_",id)]]
          if(is.null(input[[paste0(names[i],"_",id)]]) | input[[paste0(names[i],"_",id)]] == ""){
            showNotification("Names must be provided for each group", type = 'warning')
            return()
          }
        }
      }else{
        if(all(rvalues$tmp_annotations=='unlabled')){
          showNotification("No annotations have been added", type = 'warning')
        }else if(is.null(input$new_scheme_name) | input$new_scheme_name==""){
          showNotification("Please enter an name for the annotation scheme", type = 'warning')
        }else{
          new <- rvalues$tmp_annotations
        }
      }
      for(i in 1:length(rvalues$h5ad)){
        rvalues$h5ad[[i]]$obs[input$new_scheme_name] <- new
      }
      rvalues$obs <- rvalues$h5ad[[1]]$obs_keys()
      rvalues$obs_cat <- check_if_obs_cat(rvalues$h5ad[[1]]$obs)
      showNotification("Saving...", duration = NULL, id = 'save_annot')
          for(i in 1:length(rvalues$path_to_data)){
            rvalues$h5ad[[i]]$write(filename = rvalues$path_to_data[i])
          }
      Sys.sleep(1)
      removeNotification(id = 'save_annot')
      Sys.sleep(1)
      showNotification("New annotations added!", type = 'message')
    }
  })
  
  #-- genes to search from for ncbi query --#
  output$gene_query <- renderUI({
    req(rvalues$features)
    selectInput(inputId = 'gene_query', label = 'Select a gene of interest to learn more', choices = rvalues$features, selected = NA, multiple = FALSE)
  })

  #-- get gene summary from ncbi --#
  observeEvent(input$query_ncbi,{
    if(is.null(input$gene_query) | is.na(input$gene_query) | identical(input$gene_query, "") | identical(input$gene_query, " ")){
      showNotification("A valid gene must be entered", type = "error")
      return(NULL)
    }
    id <- fromJSON(file = paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=", input$gene_query, "[sym]+AND+",input$organism,"[orgn]&retmode=json"))$esearchresult$idlist
    if(identical(id, list())){
      output$gene_summary <- renderText({"Error: No genes matching this query were found on NCBI"})
    }else{
    if(length(id)>1){
      showNotification('Caution: More that one gene matched this query. Only showing first match', type = 'warning')
      id <- id[1]
    }
    output$gene_summary <- renderText({
      results <- fromJSON(file = paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=",id,"&retmode=json"))$result[[2]]
      name <- results$name
      alias <- results$otheraliases
      summary <- results$summary
      return(paste0("[id=",id,"][name=",name,"][aliases=", paste(alias, collapse = " "),"][summary=",summary,"]"))
    })
    }
  })

  ###==============// MODALITIES TAB //===================####
  shinyFileChoose(input, "h5ad_in_mod", roots = volumes, session = session)

  observeEvent(input$h5ad_in_mod, {
    path <- parseFilePaths(selection = input$h5ad_in_mod, roots = volumes)$datapath
    if(is.integer(path[1]) || identical(path, character(0)) || identical(path, character(0))) return(NULL)
    h5ad_files <- path#paste0(path,"/",list.files(path))
    assays <- sub(".h5ad","",sub(paste0(".*/"),"",h5ad_files))
    data <- list()
    ## Iterate over all assays and connect to h5ad objects
    for(i in 1:length(assays)){
      data[[i]] <- tryCatch({
        anndata$read(h5ad_files[i])
      },
      error = function(e){
        showModal(modalDialog(p(paste0("An error occured trying to connect to ", h5ad_files[i])), title = "Error connecting to h5ad file."), session = getDefaultReactiveDomain())
        return(NULL)
      })
    }
    if(is.null(data)) return(NULL)
    if(length(data) != length(assays)) return(NULL)
    if(length(unlist(lapply(data, function(x){x}))) != length(assays)) return(NULL)
    
    if(all(data[[1]]$obs_names$values == rvalues$h5ad[[1]]$obs_names$values) == FALSE){
        showModal(modalDialog(p(paste0("Cell mismatch error. Not all cell IDs from ", h5ad_files[i], " match the cell IDs from ", rvalues$path_to_data, ".")), title = "Error connecting to h5ad file."), session = getDefaultReactiveDomain())
        return(NULL)
    }

    names(data) <- assays
    
    if(is.null(data[[1]]$raw$var)){
      rvalues_mod$features <- rownames(data[[1]]$var)
    }else{
      test_gene_name <- rownames(data[[1]]$var)[1]
      if(test_gene_name %in% rownames(data[[1]]$raw$var)){ # check if rownames are numbers or gene names
        rvalues_mod$features <- rownames(data[[1]]$raw$var)
      }else if("features" %in% colnames(data[[1]]$raw$var)){ ## Check if there is a column named features in raw
        rvalues_mod$features <- data[[1]]$raw$var$features
      }else if(test_gene_name %in% data[[1]]$raw$var[,1]){ # otherwise, check if the first column contains rownames
        rvalues_mod$features <- data[[1]]$raw$var[,1]
      }
    }
    rvalues_mod$obs <- data[[1]]$obs_keys()
    ## Determine type of annotation and create a layer to annotate for easy usage later on
    rvalues_mod$obs_cat <- check_if_obs_cat(obs_df = data[[1]]$obs) ## Function to check if an observation is categorical or numeric
    reductions <- data[[1]]$obsm$as_dict()
    reduction_keys <- data[[1]]$obsm_keys()
    r_names <- rownames(data[[1]]$obs)
    for(i in 1:length(reductions)){
      reductions[[i]] <- as.data.frame(reductions[[i]])
      colnames(reductions[[i]]) <- paste0(reduction_keys[i], "_", 1:ncol(reductions[[i]]))
      rownames(reductions[[i]]) <- r_names
    }
    names(reductions) <- reduction_keys
    rvalues_mod$reductions <- reductions
    rvalues_mod$cell_ids <- rownames(data[[1]]$obs)
    rvalues_mod$h5ad <- data
    rvalues_mod$path_to_data <- h5ad_files
    init <<- 0
  })

  #-- select input for Assay --#
  output$assay_mod <- renderUI({
    req(rvalues_mod$h5ad)
    selectInput('assay_mod', "Select Assay", 
      choices = names(rvalues_mod$h5ad), 
      selected = ifelse(any(names(rvalues_mod$h5ad)=="RNA"),
        yes = "RNA",
        no = names(rvalues_mod$h5ad)[1])
      )
  })

  #-- Display name of data --#
  output$data_used_mod <- renderText({
    if(!is.null(input$assay_mod) & !is.na(input$h5ad_in_mod[[1]][2])){
      path <- sub(".*\\/", "", parseFilePaths(selection = input$h5ad_in_mod, roots = volumes)$datapath)
      return(paste0("Chosen data: ", path))
    }else{
      return(paste0(""))
    }
  })

  output$grouping_mod <- renderUI({
    req(input$assay_mod, rvalues_mod$obs)
    assay <- input$assay_mod
    if(any(grepl("seurat_clusters", rvalues_mod$obs, ignore.case = TRUE))){
      sel <- rvalues_mod$obs[grep("seurat_clusters", rvalues_mod$obs, , ignore.case = TRUE)]
    }else if(any(grepl(paste0(tolower(assay),"_clusters"), rvalues_mod$obs, ignore.case = TRUE))){
      sel <- rvalues_mod$obs[grep(paste0(tolower(assay),"_clusters"), rvalues_mod$obs, , ignore.case = TRUE)]
    }else if(any(grepl("louvain", rvalues_mod$obs, ignore.case = TRUE))){
      sel <- rvalues_mod$obs[grep("louvain", rvalues_mod$obs, ignore.case = TRUE)]
    }else{
      sel <- rvalues_mod$obs[1]
    }
    selectInput(inputId = 'grouping_mod', label = 'Group By', choices = rvalues_mod$obs, selected = sel, multiple = FALSE)
  })

  output$reduction_mod <- renderUI({
    req(input$assay_mod)
    assay <- input$assay_mod
    options <- names(rvalues_mod$reductions)
    sel <- if(any(grepl("umap", options, ignore.case = TRUE))){
      options[grepl("umap", options, ignore.case = TRUE)][1]
    }else if(any(grepl("tsne", options, ignore.case = TRUE))){
      options[grepl("tsne", options, ignore.case = TRUE)][1]
    }else{
      options[1]
    }
    selectInput(inputId = 'reduction_mod', 'Choose Clustering Method', choices = as.list(options), selected = sel)
  })

  output$featureplot_mod_feature_select <- renderUI({
    req(input$assay_mod)
    assay <- input$assay_mod
    selectInput(inputId = 'featureplot_mod_feature_select', 
                label = 'Select a Feature to Visualize on the Feature Plot', 
                choices = rvalues_mod$features, 
                selected = rvalues_mod$features[1], 
                multiple = ifelse(input$nebulosa_on == "yes", TRUE, FALSE))
  })

  output$ridgeplot_mod_feature_select <- renderUI({
    req(input$assay_mod)
    assay <- input$assay_mod
    selectInput(inputId = 'ridgeplot_mod_feature_select', 
                label = 'Select Features to Visualize on the Ridge Plot', 
                choices = rvalues_mod$features, 
                selected = rvalues_mod$features[1], 
                multiple = TRUE)
  })

  #-- dimensional reduction plot coloured by cell groups --#
  output$dimplot_mod <- renderPlotly({
    req(input$grouping_mod, input$reduction_mod)
    group.by <- list(rvalues_mod$h5ad[[1]]$obs[input$grouping_mod][,,drop=TRUE])
    names(group.by) <- input$grouping_mod
    names(group.by[[1]]) <- rvalues_mod$cell_ids
    cat <- rvalues_mod$obs_cat[which(rvalues_mod$obs == input$grouping_mod)]
    if(cat){
      dimPlotlyOutput(assay.in = input$assay_mod, 
                    reduc.in = rvalues_mod$reductions[[input$reduction_mod]], 
                    group.by = group.by, 
                    annot_panel = "", 
                    low.res = 'yes')
    }else{
      feature.in <- rvalues_mod$h5ad[[1]]$obs[input$grouping_mod][,,drop=FALSE]
      colnames(feature.in) <- input$grouping_mod
      rownames(feature.in) <- rownames(rvalues_mod$reductions[[input$reduction_mod]])
      featurePlotlyOutput(assay.in = input$assay_mod,
                          reduc.in = rvalues_mod$reductions[[input$reduction_mod]],
                          group.by = group.by,
                          feature.in = feature.in,
                          low.res = 'yes')
    }
  })

  #-- dimensional reduction plot coloured by feature expression --#
  output$featureplot_mod <- renderPlotly({
    req(input$grouping_mod, input$reduction_mod, input$featureplot_mod_feature_select)
    group.by <- list(rvalues_mod$h5ad[[1]]$obs[input$grouping_mod][,,drop=TRUE])
    names(group.by) <- input$grouping_mod
    names(group.by[[1]]) <- rvalues_mod$cell_ids
    if(rvalues_mod$raw_dtype == "NULL" | rvalues$raw_dtype == "counts"){
      feature.in <- as.data.frame(rvalues_mod$h5ad[[1]]$X[,match(input$featureplot_mod_feature_select, rvalues_mod$features)])
    }else if(rvalues_mod$raw_dtype == "normalized"){
      feature.in <- as.data.frame(rvalues_mod$h5ad[[1]]$raw$X[,match(input$featureplot_mod_feature_select, rvalues_mod$features)])
    }
    colnames(feature.in) <- input$featureplot_mod_feature_select
    rownames(feature.in) <- rownames(rvalues_mod$reductions[[input$reduction_mod]])
    if(input$nebulosa_mod_on == 'no'){
      featurePlotlyOutput(assay.in = input$assay_mod,
                          reduc.in = rvalues_mod$reductions[[input$reduction_mod]],
                          group.by = group.by,
                          feature.in = feature.in,
                          low.res = 'yes')
    }else{
      featurePlotlyOutput_nebulosa(assay.in = input$assay_mod,
                                reduc.in = rvalues_mod$reductions[[input$reduction_mod]],
                                group.by = as.data.frame(group.by),
                                feature.in = feature.in,
                                low.res = 'yes')
    }
  })

  ## output variable to hold type of user selected grouping for main panel
  output$grouping_mod_type <- reactive({
    req(input$grouping_mod)
    cat <- rvalues_mod$obs_cat[which(rvalues_mod$obs == input$grouping_mod)]
    if(cat){"yes"}else{"no"}
  })
  outputOptions(output, "grouping_mod_type", suspendWhenHidden = FALSE)

  #-- dimensional reduction plot coloured by cell groups --#
  output$ridgeplot_mod <- renderPlot({
    req(input$ridgeplot_mod_feature_select, input$grouping_mod)
    if(is.null(rvalues_mod$h5ad[[1]]$raw)){
      data.features <- as.data.frame(rvalues_mod$h5ad[[1]]$X[,match(input$ridgeplot_mod_feature_select, rvalues_mod$features)])
    }else{
      data.features <- as.data.frame(rvalues_mod$h5ad[[1]]$raw$X[,match(input$ridgeplot_mod_feature_select, rvalues_mod$features)])
    }
    colnames(data.features) <- input$ridgeplot_mod_feature_select
    rownames(data.features) <- rownames(rvalues_mod$reductions)
    data.features$id <- as.character(rvalues_mod$h5ad[[1]]$obs[input$grouping_mod][,,drop=TRUE])

    data.features <- reshape2::melt(data.features)
    colnames(data.features) <- c("ident", "feature", "expression")

    ggRidgePlot(data.features)
  })
  
  #### CRISPR feature 1 parameters
  output$crispr_feature_1 <- renderUI({ ## Feature 1 selected for CRISPR
    req(input$assay_mod)
    selectInput(inputId = 'crispr_feature_1_sel', 
                label = 'Select feature #1 you want to compare!', 
                choices = rvalues_mod$features, 
                selected = rvalues_mod$features[1], 
                multiple = FALSE)
  })
  
  output$crispr_feature_1_slider <- renderUI({ ## Feature 1 selected for CRISPR
    req(input$assay_mod)
    req(input$crispr_feature_1_sel)
    if(rvalues_mod$raw_dtype == "NULL" | rvalues$raw_dtype == "counts"){
      expr_range <- rvalues_mod$h5ad[[1]]$X[,match(input$crispr_feature_1_sel, rvalues_mod$features)]
    }else if(rvalues_mod$raw_dtype == "normalized"){
      expr_range <- rvalues_mod$h5ad[[1]]$raw$X[,match(input$crispr_feature_1_sel, rvalues_mod$features)]
    }
    min_val <- round(min(expr_range),2)
    max_val <- round(max(expr_range),2)
    sel_value <- median((expr_range))
    sliderInput(inputId = 'crispr_feature_1_slider_val', 
                label = 'Set your expression cutoff for feature 1!', 
                min = min_val,
                max = max_val,
                value = sel_value)
  })
  
  output$crispr_feature_1_dist <- renderPlot({
    req(input$assay_mod)
    req(input$crispr_feature_1_sel)
    req(input$crispr_feature_1_slider_val)
    
    if(is.null(rvalues_mod$h5ad[[1]]$raw)){
      crispr_exp_feature_1_df <-  data.frame("exp" =rvalues_mod$h5ad[[1]]$X[,match(input$crispr_feature_1_sel, rvalues_mod$features)],
                                             "feature" = input$crispr_feature_1_sel)
    }else{
      crispr_exp_feature_1_df <-  data.frame("exp" =rvalues_mod$h5ad[[1]]$raw$X[,match(input$crispr_feature_1_sel, rvalues_mod$features)],
                                             "feature" = input$crispr_feature_1_sel)
    }

    ggplot(crispr_exp_feature_1_df,aes(exp,fill = feature)) +
      geom_density(fill = "#4682B4") +
      theme_cowplot() +
      geom_vline(xintercept = input$crispr_feature_1_slider_val,
                 size = 1.5,
                 color = "black",
                 linetype = 2) +
      labs(title = input$crispr_feature_1_sel)
  })
  
  crispr_feature_1_cells <- reactive({
    req(input$assay_mod)
    req(input$crispr_feature_1_sel)
    req(input$crispr_feature_1_slider_val)

    if(rvalues_mod$raw_dtype == "NULL" | rvalues$raw_dtype == "counts"){
      crispr_exp_feature_1_df <-  data.frame("exp" =rvalues_mod$h5ad[[1]]$X[,match(input$crispr_feature_1_sel, rvalues_mod$features)],
                                             "feature" = input$crispr_feature_1_sel)
    }else if(rvalues_mod$raw_dtype == "normalized"){
      crispr_exp_feature_1_df <-  data.frame("exp" =rvalues_mod$h5ad[[1]]$raw$X[,match(input$crispr_feature_1_sel, rvalues_mod$features)],
                                             "feature" = input$crispr_feature_1_sel)
    }

    crispr_exp_feature_1_df <- crispr_exp_feature_1_df %>%
      mutate("pass_exp_thr" = if_else(exp >= input$crispr_feature_1_slider_val,"yes","no"),
             "index" = rownames(crispr_exp_feature_1_df))
    crispr_exp_feature_1_df$cell_id <- rvalues_mod$cell_ids
    crispr_exp_feature_1_df
  })
  
  output$crispr_feature_1_cells_print <- renderText({
    req(crispr_feature_1_cells())
    cells_pass <- subset(crispr_feature_1_cells(),pass_exp_thr == "yes")
    res_string <- paste("There are ",nrow(cells_pass)," cells above your selected threshold for feature:",
                        input$crispr_feature_1_sel,sep="")
    res_string
  })

  
  #### CRISPR feature 2 parameters
  output$crispr_feature_2 <- renderUI({ ## Feature 1 selected for CRISPR
    req(input$assay_mod)
    req(input$crispr_feature_1_sel)
    
    feature_options <- setdiff(rvalues_mod$features,input$crispr_feature_1_sel)

    selectInput(inputId = 'crispr_feature_2_sel', 
                label = 'Select feature #2 you want to compare!', 
                choices = feature_options, 
                selected = feature_options[1], 
                multiple = FALSE)
  })
  
  output$crispr_feature_2_slider <- renderUI({ ## Feature 1 selected for CRISPR
    req(input$assay_mod)
    req(input$crispr_feature_2_sel)
    if(rvalues_mod$raw_dtype == "NULL" | rvalues$raw_dtype == "counts"){
      expr_range <- rvalues_mod$h5ad[[1]]$X[,match(input$crispr_feature_2_sel, rvalues_mod$features)]
    }else if(rvalues_mod$raw_dtype == "normalized"){
      expr_range <- rvalues_mod$h5ad[[1]]$raw$X[,match(input$crispr_feature_2_sel, rvalues_mod$features)]
    }
    min_val <- round(min(expr_range),2)
    max_val <- round(max(expr_range),2)
    sel_value <- median(expr_range)
    sliderInput(inputId = 'crispr_feature_2_slider_val', 
                label = 'Set your expression cutoff for feature 2!', 
                min = min_val,
                max = max_val,
                value = sel_value)
  })
  
  output$crispr_feature_2_dist <- renderPlot({
    req(input$assay_mod)
    req(input$crispr_feature_2_sel)
    req(input$crispr_feature_2_slider_val)
    
    if(rvalues_mod$raw_dtype == "NULL" | rvalues$raw_dtype == "counts"){
      crispr_exp_feature_2_df <-  data.frame("exp" =rvalues_mod$h5ad[[1]]$X[,match(input$crispr_feature_2_sel, 
                                                                                       rvalues_mod$features)],
                                             "feature" = input$crispr_feature_2_sel)
    }else if(rvalues_mod$raw_dtype == "normalized"){
      crispr_exp_feature_2_df <-  data.frame("exp" =rvalues_mod$h5ad[[1]]$raw$X[,match(input$crispr_feature_2_sel, 
                                                                                       rvalues_mod$features)],
                                             "feature" = input$crispr_feature_2_sel)
    }
  
    ggplot(crispr_exp_feature_2_df,aes(exp,fill = feature)) +
      geom_density(fill = "#B47846") +
      theme_cowplot() +
      geom_vline(xintercept = input$crispr_feature_2_slider_val,
                 size = 1.5,
                 color = "black",
                 linetype = 2) +
      labs(title = input$crispr_feature_2_sel)
    
  })
  
  crispr_feature_2_cells <- reactive({
    req(input$assay_mod)
    req(input$crispr_feature_2_sel)
    req(input$crispr_feature_2_slider_val)
    
    if(rvalues_mod$raw_dtype == "NULL" | rvalues$raw_dtype == "counts"){
      crispr_exp_feature_2_df <-  data.frame("exp" =rvalues_mod$h5ad[[1]]$X[,match(input$crispr_feature_2_sel, rvalues_mod$features)],
                                             "feature" = input$crispr_feature_2_sel)
    }else if(rvalues_mod$raw_dtype == "normalized"){
      crispr_exp_feature_2_df <-  data.frame("exp" =rvalues_mod$h5ad[[1]]$raw$X[,match(input$crispr_feature_2_sel, rvalues_mod$features)],
                                             "feature" = input$crispr_feature_2_sel)
    }

    crispr_exp_feature_2_df <- crispr_exp_feature_2_df %>%
      mutate("pass_exp_thr" = if_else(exp >= input$crispr_feature_2_slider_val,"yes","no"),
             "index" = rownames(crispr_exp_feature_2_df))
    crispr_exp_feature_2_df$cell_id <- rvalues_mod$cell_ids
    crispr_exp_feature_2_df
  })
  
  output$crispr_feature_2_cells_print <- renderText({
    req(crispr_feature_2_cells())
    cells_pass <- subset(crispr_feature_2_cells(),pass_exp_thr == "yes")
    res_string <- paste("There are ",nrow(cells_pass)," cells above your selected threshold for feature:",
                        input$crispr_feature_2_sel,sep="")
    res_string
  })
  
  #### Action button for calculating CRISPR DE
  observeEvent(input$crispr_de_analysis,{
    shinyjs::showElement(id= "crispr_res")
    
    crispr_feature_1_cells_rna <- isolate({
      #req(crispr_feature_1_cells())
      pos_cells <- subset(crispr_feature_1_cells(),pass_exp_thr == "yes")
      if(is.null(rvalues_mod$h5ad[[1]]$raw)){
        hvg_features <- rownames(rvalues$h5ad[[1]]$var)
        crispr_feature_1_rna_mat <- as.data.frame(rvalues$h5ad[[1]]$X[match(pos_cells$cell_id,rvalues$cell_ids),hvg_features])
      }else{
        hvg_features <- match(rownames(rvalues$h5ad[[1]]$var),rownames(rvalues$h5ad[[1]]$raw$var))
        crispr_feature_1_rna_mat <- as.data.frame(rvalues$h5ad[[1]]$raw$X[match(pos_cells$cell_id,rvalues$cell_ids),hvg_features])
      }
      crispr_feature_1_rna_mat
    })
    
    crispr_feature_2_cells_rna <- isolate({
      #req(crispr_feature_2_cells())
      pos_cells <- subset(crispr_feature_2_cells(),pass_exp_thr == "yes")
      if(rvalues_mod$raw_dtype == "NULL" | rvalues$raw_dtype == "counts"){
        hvg_features <- rownames(rvalues$h5ad[[1]]$var)
        crispr_feature_2_rna_mat <- as.data.frame(rvalues$h5ad[[1]]$X[match(pos_cells$cell_id,rvalues$cell_ids),hvg_features])
      }else if(rvalues_mod$raw_dtype == "normalized"){
        hvg_features <- match(rownames(rvalues$h5ad[[1]]$var),rownames(rvalues$h5ad[[1]]$raw$var))
        crispr_feature_2_rna_mat <- as.data.frame(rvalues$h5ad[[1]]$raw$X[match(pos_cells$cell_id,rvalues$cell_ids),hvg_features])
      }

      crispr_feature_2_rna_mat
    })
    
    merged_crispr_feature_avg_exp <- reactive({
      
      crispr_feature_1_rna_mat_means <- colMeans(as.matrix(crispr_feature_1_cells_rna))
      crispr_feature_2_rna_mat_means <- colMeans(as.matrix(crispr_feature_2_cells_rna))
      gene_names <- rownames(rvalues$h5ad[[1]]$var)
      
      merged_colMeans <- data.frame("feature_1_avg_exp"= crispr_feature_1_rna_mat_means,
                                    "feature_2_avg_exp" = crispr_feature_2_rna_mat_means,
                                    "gene" = gene_names)
      
      merged_colMeans <- merged_colMeans %>%
        mutate("diff_features_avg_exp" = feature_1_avg_exp - feature_2_avg_exp) %>%
        mutate_if(is.numeric, round,2)
      
      merged_colMeans
    })
      
    output$crispr_avg_gene_exp <- renderPlotly({
      req(merged_crispr_feature_avg_exp())
      avg_exp_plot_plotly <- plot_ly(data = merged_crispr_feature_avg_exp(),
                                     x = ~feature_1_avg_exp,
                                     y = ~feature_2_avg_exp,
                                     text = ~paste("Gene: ", gene,sep=" "),
                                     type = 'scatter',
                                     mode = 'markers',
                                     marker = list(size = 6,
                                                   color = 'black'))
      avg_exp_plot_plotly
    })

    selected <- reactive(getReactableState("crispr_avg_gene_exp_tbl", "selected"))
    
    output$crispr_avg_gene_exp_tbl <- renderReactable({
      req(merged_crispr_feature_avg_exp())
    
      reactable(merged_crispr_feature_avg_exp(),
                selection = 'single' , onClick = 'select', searchable = TRUE,)
    })
    
    output$crispr_gene_vlnplot <- renderPlot({
      req(merged_crispr_feature_avg_exp())
      req(selected())
      
      gene_selected <- merged_crispr_feature_avg_exp()[selected(),]$gene
      gene_selected_index <- match(gene_selected,rownames(rvalues$h5ad[[1]]$var))
      
      feature_1_gene_exp <- crispr_feature_1_cells_rna
      feature_1_gene_exp <- feature_1_gene_exp[,gene_selected_index]
      feature_1_gene_exp_df <- data.frame("exp" = feature_1_gene_exp,
                                            "group" = isolate(input$crispr_feature_1_sel))
      feature_2_gene_exp <- crispr_feature_2_cells_rna
      feature_2_gene_exp <- feature_2_gene_exp[,gene_selected_index]
      feature_2_gene_exp_df <- data.frame("exp" = feature_2_gene_exp,
                                          "group" =  isolate(input$crispr_feature_2_sel))
      merged_feature_gene_exp <- rbind(feature_1_gene_exp_df,feature_2_gene_exp_df)

      ggplot(merged_feature_gene_exp,aes(group,exp,fill = group)) +
        geom_violin() +
        stat_summary(fun=mean, geom="point", size=5, color = "black") +
        scale_fill_manual(values = c("#4682B4","#B47846")) +
        theme_cowplot() +
        theme(legend.position = "none") +
        labs(x = "Features",
             y = "Gene expression",
             title = gene_selected)
    })
    
  })

  ###==============// CUSTOM META DATA TAB //==============####

  #-- select input for Assay --#
  output$assay_3 <- renderUI({
    req(rvalues$h5ad)
    selectInput('assay_3', "Select Assay", 
      choices = names(rvalues$h5ad), 
      selected = ifelse(any(names(rvalues$h5ad)=="RNA"),
        yes = "RNA",
        no = names(rvalues$h5ad)[1])
      )  
  })

  #-- select meta data to combine --#
  output$meta_group_checkbox <- renderUI({
    req(input$assay_3, rvalues$h5ad)
    checkboxGroupInput(inputId = 'meta_group_checkbox', label = 'Meta Data', choices = rvalues$obs)
  })

  #-- display the order of the selected meta data --#
  output$example_meta_group_name <- renderText({
    req(rvalues$h5ad)
    if(is.null(rvalues$order) || is.na(rvalues$order) || length(rvalues$order) == 0){
      return(NULL)
    }
    paste(rvalues$order,collapse=" + ")
  })

  #-- display a sample of what the combined meta data will resemble --#
  output$example_meta_group <- renderTable({
    req(input$assay_3)
    if(is.null(rvalues$order) || is.na(rvalues$order) || length(rvalues$order) == 0){
      return(NULL)
    }
    ex <- rvalues$h5ad[[1]]$obs[rvalues$order]
    ex_df <- as.data.frame(apply(X = ex, MARGIN = 1, FUN = paste, collapse = '_'), stringsAsFactors = FALSE)[sample(x = 1:nrow(ex),size = 20,replace = FALSE),,drop=F]
    colnames(ex_df) <- 'Examples of chosen meta data combination'
    return(ex_df)
  }, hover = TRUE, width = '100%',align = 'c')

  #-- update the 'order' global variable whenever the checkbox group is selected/deselected --#
  observeEvent(input$meta_group_checkbox, ignoreNULL = F,{
    req(rvalues$h5ad)
    if(!is.null(input$meta_group_checkbox)){
      new <- input$meta_group_checkbox
      if(length(new)>length(rvalues$order) | is.null(rvalues$order)){
        new <- new[!new%in%rvalues$order]
        rvalues$order <- c(rvalues$order, new)
      }else if(identical(new, character(0))){
        rvalues$order <- NULL
      }else{
        rvalues$order <- rvalues$order[-which(!rvalues$order%in%new)]
      }
    }else{
    rvalues$order <- NULL
    }
  })

  #-- add custom metadata grouping to loom file --#
  observeEvent(input$add_new_meta_group,{
    if(is.null(input$meta_group_checkbox) || length(input$meta_group_checkbox)<2){
      showNotification('At least two meta data groups must be selected!', type = 'error')
    }else if(is.null(input$new_meta_group) || input$new_meta_group == ''){
      showNotification('A name for the meta data must be provided!', type = 'error')
    }else if(input$new_meta_group %in% rvalues$obs){
      showNotification('The chosen meta data name already exists! Please select a unique name.', type = 'error')
    }else{
      ex <- rvalues$h5ad[[1]]$obs[rvalues$order]
      ex_list <- list(apply(X = ex, MARGIN = 1, FUN = paste, collapse = '_'))
      showModal(modalDialog(p(paste0("Adding ", input$new_meta_group, " to meta data...")), title = "This window will close after completion"), session = getDefaultReactiveDomain())
      Sys.sleep(2)
      for(i in 1:length(rvalues$h5ad)){
        rvalues$h5ad[[i]]$obs[input$new_meta_group] <- unlist(ex_list)
      }
      Sys.sleep(2)
      removeModal(session = getDefaultReactiveDomain())
      showNotification("Saving...", duration = NULL, id = 'save_meta')
      rvalues$obs <- rvalues$h5ad[[1]]$obs_keys()
      rvalues$obs_cat <- check_if_obs_cat(rvalues$h5ad[[1]]$obs)
      for(i in 1:length(rvalues$path_to_data)){
        rvalues$h5ad[[i]]$write(filename = rvalues$path_to_data[i])
      }
      removeNotification(id = 'save_meta')
    }
  })

  ###==============// FILE CONVERSION TAB //==============####
  
  shinyFileChoose(input, "object_file", roots = volumes, session = session)
  shinyDirChoose(input, "new_directory", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  output$chosen_object_file <- renderText(parseFilePaths(selection = input$object_file, roots = volumes)$datapath)
  
  output$chosen_new_directory <- renderText(parseDirPath(roots = volumes, selection = input$new_directory))
  
  observeEvent(input$object_file, ignoreInit = T, ignoreNULL = T, {
    file_path <- parseFilePaths(selection = input$object_file, roots = volumes)$datapath
    if(identical(file_path,character(0))){
      return(NULL)
    }
    file_size <- fileSize(file_path, units = 'GB')
    if(file_size>1){
      output$notes <- renderText(paste0('Warning: The selected file is larger than 1 GB (',round(file_size,2),' GB) and may take some time to load depending on your computer (if running locally) or the server (if running on the cloud)'))
      }else{
        output$notes <- NULL
      }
      })
  
  
  observeEvent(input$sl_convert, ignoreInit = T, {
      if(is.integer(input$object_file[1])){
        showNotification('A Seurat or Scanpy object must be selected', type = 'error')
        return(NULL)
      }else if(is.integer(input$new_directory[1])){
        showNotification('A directory to save the converted file to must be selected', type = 'error')
        return(NULL)
      }
      dir_path <- parseDirPath(roots = volumes, selection = input$new_directory)

      if(length(list.files(dir_path))>0){
        showNotification('The selected directory already contains files. Please used a new/empty directory.', type = 'error')
        return(NULL)
      }

      file_path <- parseFilePaths(selection = input$object_file, roots = volumes)$datapath

      if(!grepl('\\.rds$', file_path, ignore.case = T)){
        if(!grepl('\\.h5|\\.h5ad', file_path, ignore.case = T)){
          showNotification('The selected file must be of type .rds or .h5/.h5ad', type = 'error')
          return(NULL)
        }
      }

      showModal(modalDialog(p("Converting Seurat Object to loom file(s). Please Wait..."), title = "This window will close after conversion is complete"), session = getDefaultReactiveDomain())
      error <- seuratToLoom(obj = file_path, dir = dir_path)
      removeModal(session = getDefaultReactiveDomain())
      if(error != 1){
        showNotification('Error in file conversion!', type = 'error', closeButton = TRUE, duration = 100)
      }else{
        showNotification('Conversion Complete!', type = 'message', closeButton = F, duration = 15)
      }
    })
  
  shinyFileChoose(input, "update_file", roots = volumes, session = session)
  shinyDirChoose(input, "new_directory_2", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  output$chosen_update_file <- renderText(parseFilePaths(selection = input$update_file, roots = volumes)$datapath)
  
  output$chosen_new_directory_2 <- renderText(parseDirPath(roots = volumes, selection = input$new_directory_2))
  
  observeEvent(input$ls_convert, ignoreInit = T, {
      if(is.integer(input$update_file[1])){
        showNotification('A seurat object must be selected', type = 'error')
        return(NULL)
      }else if(is.integer(input$new_directory_2)){
        showNotification('A directory to save the converted file to must be selected', type = 'error')
        return(NULL)
      }else if(input$new_file_name == "" & input$overwrite == "no"){
        showNotification('A file name for the updated objet must be provided', type = 'error')
        return(NULL)
      }
      dir_path <- parseDirPath(roots = volumes, selection = input$new_directory_2)
      path_to_seurat <- parseFilePaths(selection = input$update_file, roots = volumes)$datapath

      if(input$overwrite == "no"){
        file <- input$new_file_name
        if(grepl(".rds$", file, ignore.case = T) == F){
          file <- paste0(file,".rds")
        }
      }else{
        file <- NULL
      }

      if(!grepl('\\.rds$', path_to_seurat, ignore.case = T)){
        showNotification('The selected file must be of type .rds', type = 'error')
        return(NULL)
      }

      showModal(modalDialog(p("Converting Loom to Seurat. Please Wait..."), title = "This window will close after conversion is complete"), session = getDefaultReactiveDomain())
      error <- loomToSeurat(obj = path_to_seurat, loom = loom_data(), dir = dir_path  , file = file)
      removeModal(session = getDefaultReactiveDomain())
      if(error != 1){
        showNotification('Error in file conversion!', type = 'error', closeButton = TRUE, duration = 100)
      }else{
        showNotification('Conversion Complete!', type = 'message', closeButton = F, duration = 15)
      }
    })

  shinyFileChoose(input, "in_file_to_h5ad", roots = volumes, session = session)
  shinyDirChoose(input, "out_dir_to_h5ad", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyFileChoose(input, "legacy_file_to_h5ad", roots = volumes, session = session)

  output$chosen_in_file_to_h5ad <- renderText(parseFilePaths(selection = input$in_file_to_h5ad, roots = volumes)$datapath)
  output$chosen_out_file_to_h5ad <- renderText(input$out_file_to_h5ad)
  output$chosen_out_dir_to_h5ad <- renderText(parseDirPath(roots = volumes, selection = input$out_dir_to_h5ad))
  output$chosen_legacy_file_to_h5ad <- renderText(parseFilePaths(selection = input$legacy_file_to_h5ad, roots = volumes)$datapath)

  observeEvent(input$to_h5ad, ignoreInit = T, {

    if(is.integer(input$in_file_to_h5ad[1])){
      showNotification('A file to convert must be selected', type = 'error')
      return(NULL)
    }else if(input$out_file_to_h5ad == ""){
      showNotification('An output h5ad file name must be specified', type = 'error')
      return(NULL)
    }else if(is.integer(input$out_dir_to_h5ad[1])){
      showNotification('A directory to save the converted file to must be selected', type = 'error')
      return(NULL)
    }else if(is.integer(input$legacy_file_to_h5ad[1]) & input$is_legacy == "yes"){
      showNotification('The associated seurat object to the legacy loom file must be selected', type = 'error')
      return(NULL)
    }
    
    dir_path <- parseDirPath(roots = volumes, selection = input$out_dir_to_h5ad)
    path_to_in_file <- parseFilePaths(selection = input$in_file_to_h5ad, roots = volumes)$datapath

    if(input$is_legacy == "yes"){
      path_to_legacy <- parseFilePaths(selection = input$legacy_file_to_h5ad, roots = volumes)$datapath
    }else{
      path_to_legacy <- NULL
    }

    out_file <- input$out_file_to_h5ad
    if(grepl("\\.h5ad$", out_file) == FALSE){
      out_file <- paste0(out_file, ".h5ad")
    }
    out_path <- file.path(dir_path, out_file)
    if(out_file %in% list.files(dir_path)){
      showModal(modalDialog(p(paste0("Error: ", out_file, " already exisits in ", dir_path, ". Please choose a unique file.")), title = "File already exists."), session = getDefaultReactiveDomain())
      return(NULL)
    }

    showModal(modalDialog(p("Converting to h5ad. Please Wait..."), title = "This window will close after conversion is complete"), session = getDefaultReactiveDomain())
    error <- scap_to_h5ad(in_file = path_to_in_file, out_path = out_path, old_file = path_to_legacy)
    removeModal(session = getDefaultReactiveDomain())
    if(error == 0){
      showModal(modalDialog(p("Expecting Rds or loom file"), title = "Error"), session = getDefaultReactiveDomain())
    }else if(error == -1){
      showModal(modalDialog(p("The loom file you are trying to convert is in a legacy format. Please specify the associated Seurat Object and try again."), title = "Error"), session = getDefaultReactiveDomain())
    }else if(error == -2){
      showModal(modalDialog(p("Legacy conversion error. Could not proceed."), title = "Error"), session = getDefaultReactiveDomain())
    }else if(error == -3){
      showModal(modalDialog(p("The associated Seurat object you selected does not appear to be a Seurat object."), title = "Error"), session = getDefaultReactiveDomain())
    }else if(error == -4){
      showModal(modalDialog(p(paste0("Expecting an object of calss Seurat or loom. Instead recieved an object of class: ", class(obj)[1])), title = "Error"), session = getDefaultReactiveDomain())
    }else if(error == 1){
      showNotification('Conversion Complete!', type = 'message', closeButton = FALSE, duration = 15)
    }else{
      showNotification('Error in file conversion!', type = 'error', closeButton = TRUE, duration = 100)
    }
  })

#   ####
#   #-- Functions for annotation comparison--#
#   ####
  #-- get the first list of potential annotations from loom--#
  output$comp_anno_list1 <- renderUI({
    req(input$assay_1)
    req(rvalues$obs)
    req(rvalues$obs_cat)

    ## get annotation options from rvalues
    annotation_options <- rvalues$obs[rvalues$obs_cat]
    
    #group.by <- list(rvalues$h5ad[[1]]$obs[input$grouping_1][,,drop=TRUE])
    selectInput(
      inputId = 'comp_anno_1',
      label = 'Select the annotation you want to compare!',
      choices = annotation_options,
      multiple = FALSE)
    })

  #-- second list of annotations to compare against--#
  output$comp_anno_list2 <- renderUI({
    req(input$assay_1,input$comp_anno_1)

    annotation_options <- rvalues$obs[rvalues$obs_cat]
    annotation_options <- annotation_options[!grepl(input$comp_anno_1,annotation_options)]

    selectInput(
      inputId = 'comp_anno_2',
      label = 'Select the annotation to compare against!',
      choices = annotation_options,
      multiple = FALSE)  
    })


  ## Function to create data table for sankey diagram
  sankey_comp <- eventReactive(input$compare_annos, {
    req(input$comp_anno_1, input$comp_anno_2)

    annos_to_compare <- rvalues$h5ad[[1]]$obs[c(input$comp_anno_1, input$comp_anno_2)]

    annos_to_compare_stats <- annos_to_compare %>%
      group_by(get(input$comp_anno_1),get(input$comp_anno_2)) %>%
      tally() %>%
      ungroup()

    colnames(annos_to_compare_stats) <- c("anno1","anno2","n")
    
    annos_to_compare_stats <- as.data.frame(annos_to_compare_stats)
    rownames(annos_to_compare_stats) <- 1:nrow(annos_to_compare_stats)

    # annos_to_compare_stats$anno1 <- as.character(annos_to_compare_stats$anno1)
    # annos_to_compare_stats$anno2 <- as.character(annos_to_compare_stats$anno2)

    annos_to_compare_stats <- annos_to_compare_stats %>%
      mutate("anno1" = paste(input$comp_anno_1,anno1,sep=": ")) %>%
      mutate("anno2" = paste(input$comp_anno_2,anno2,sep=": "))

    joined_annos_1 <- unique(annos_to_compare_stats$anno1) 
    joined_annos_2 <- unique(annos_to_compare_stats$anno2)

    annos_to_compare_stats$IDsource= match(annos_to_compare_stats$anno1, joined_annos_1) - 1
    annos_to_compare_stats$IDtarget= (match(annos_to_compare_stats$anno2, joined_annos_2))  - 1
    annos_to_compare_stats <- annos_to_compare_stats %>%
      mutate("IDtarget" = IDtarget + length(joined_annos_1))

    return(annos_to_compare_stats)
    })


  #### Compare annotations elements
  ## Sankey diagram to compare two annotations
  output$sankey_diagram <- renderPlotly({
    req(sankey_comp())

    joined_annos <- c(sankey_comp()$anno1,sankey_comp()$anno2)
    joined_annos <- unique(joined_annos)

    p <- plot_ly(
      type = "sankey",
      orientation = "h",

      node = list(
        label = joined_annos,
        pad = 15,
        thickness = 40,
        line = list(
          color = "black",
          width = 0.5
          )
        ),

      link = list(
        source = sankey_comp()$IDsource,
        target = sankey_comp()$IDtarget,
        value =  sankey_comp()$n
        )
      ) %>%
      layout(font = list(
        size = 18)
        )
    
    p

    }) # sankey_diagram end

#   #### SciBet #### 
  
  ## DataTable with SciBet reference information
  output$scibet_references <- renderReactable({
    datasets_scibet_sub <- datasets_scibet %>%
      dplyr::select(Species,title,GSE,cell_types,number_of_cells,doi_link)
    
    reactable(datasets_scibet_sub,
              filterable = TRUE,
              searchable = TRUE,
              columns = list(
                Species = colDef(name = "Species", minWidth =50),
                title = colDef(name = "Title"),
                GSE = colDef(name = "GSE",  minWidth =50),
                cell_types = colDef(name = "Cell types", minWidth =200),
                number_of_cells = colDef(name = "# of cells"),
                doi_link  = colDef(name = "DOI")
                )
              )
  })
  
  ## Controls for SciBet prediction panel
  show_datasets <- reactive({
    selected_species <- input$sci_bet_species
    show_datasets <- subset(datasets_scibet,Species == selected_species)
    datasets_selected <- show_datasets$internal_id
    names(datasets_selected) <- show_datasets$title
    return(names(datasets_selected))
  })
  
  ## This output will hold the different datasets for the selected species
  output$scibet_reference_sets <- renderUI({
    req(input$sci_bet_species)
    req(show_datasets())

    selectInput(
      inputId = 'scibet_data',
      label = 'Select a dataset!',
      choices = show_datasets(),
      multiple = FALSE)
  })
  
  output$predict_cells_button <- renderUI({
    req(rvalues$h5ad)
    actionButton("predict_cells", "Predict cell types!")
  })
  
  scibet_ref_selected <- reactive({
    req(input$scibet_data)
    ref_selected <- subset(datasets_scibet,title == input$scibet_data)
    return(ref_selected)
  })
  
  #### Main function running SciBet prediction
  ## Perform predictions when user clicks button
  predictions_results <- eventReactive(input$predict_cells ,{
    ## Run prediction with the dataset selected
    req(input$assay_1)
    req(rvalues$h5ad)
    req(scibet_ref_selected())
    
    showModal(modalDialog(p(paste0("Running SciBet prediction on your data...")),
                          title = "This window will close once prediction is done! This can take a while depending on the size of your dataset!"), session = getDefaultReactiveDomain())
    #Sys.sleep(1)
    
    ## Instantiate full prediction table
    all_predictions <- data.frame()
      
    ## Get scPred model
    data_dir <- "../SciBet_reference_datasets"
    file_test_scibet <- paste(data_dir,scibet_ref_selected()$SCAP_filename,sep="/")
    model <- readr::read_csv(file_test_scibet)
    if(colnames(model)[1] == "X1" & colnames(model)[2] == "Unnamed: 0"){
      model <- model[,-1]
    }
    
    model <- pro.core(model)
    prd <- LoadModel(model)
    
    ## Read in data matrix in X chunks and predict each junk, then stitch them
    ## back together at the end
    ## For big datasets, split cell vector into chunks of ~ 2000 cells
    cell_number <- length(rvalues$cell_ids)
    sub_vector <- 1:cell_number
    n_chunks <- ceiling(cell_number/2000)
    chunk_list <- split(sub_vector, sort(sub_vector%%n_chunks))
    chunk_no <- 1
    
    for(chunk in names(chunk_list)){
      chunk_no <- chunk_no + 1
      if(is.null(rvalues$h5ad[[1]]$raw)){
        data_matrix <- rvalues$h5ad[[1]]$X[chunk_list[[chunk]],]
      }else{
        data_matrix <- rvalues$h5ad[[1]]$raw$X[chunk_list[[chunk]],]
      }
      #data_matrix <- data[[1]]$X[chunk_list[[chunk]],]
      gene.names <- rvalues$features
      colnames(data_matrix) <- gene.names
      
      ## Do prediction
      label <- prd(data_matrix)
      label_df <- data.frame("labels" = label)
      rm(data_matrix)
      
      # add prediction to full results
      all_predictions <- rbind(all_predictions,label_df)
    }
    
    #Sys.sleep(1)
    removeModal(session = getDefaultReactiveDomain())
    
    return(all_predictions)
  })
  
  ## Plot a barplot representing how many cells were assigned to which classes
  output$scibet_predictions_plot <- renderPlot({
    req(predictions_results())

    label_summary <- predictions_results() %>%
      group_by(labels) %>%
      tally() %>%
      arrange(desc(n))
    
    label_summary$labels <- factor(label_summary$labels,
                                   levels = rev(unique(label_summary$labels)))
    
    scibet_plot <- ggplot(label_summary,aes(labels,n, fill = labels)) +
      geom_bar(stat = "identity") +
      theme_cowplot() +
      theme(legend.position = "none") +
      coord_flip()
    
    return(scibet_plot)
  })
  
  output$add_predictions_button <- renderUI({
    req(rvalues$h5ad)
    req(predictions_results())
    actionButton("add_predictions_meta", "Add predictions to metadata!")
    })

  ## Add predictions to metadata
  observeEvent(input$add_predictions_meta,{
    req(predictions_results())
    req(rvalues$h5ad)
    req(input$scibet_data)
    internal_id <- subset(datasets_scibet,title == input$scibet_data)$internal_id
    pred_params_name <- paste("scibet_",internal_id,"",sep="")
    new_meta_group <- predictions_results()
    colnames(new_meta_group) <- pred_params_name

    new_meta_group_list <- as.list(new_meta_group)

    showModal(modalDialog(p(paste0("Adding predictions to meta data...")), title = "This window will close once task is completed!"), session = getDefaultReactiveDomain())
    Sys.sleep(1)
    
    rvalues$h5ad[[1]]$obs[pred_params_name] <- new_meta_group_list
    rvalues$obs <- rvalues$h5ad[[1]]$obs_keys()
    rvalues$obs_cat <- check_if_obs_cat(rvalues$h5ad[[1]]$obs)
    rvalues$h5ad[[1]]$write(filename = rvalues$path_to_data[1])
    
    Sys.sleep(1)
    removeModal(session = getDefaultReactiveDomain())
    })
  
  #### Differential expression testing elements
  output$de_annotation_list <- renderUI({
    req(input$assay_1)
    req(rvalues$obs)
    req(rvalues$obs_cat)
    
    ## get annotation options from rvalues
    annotation_options <- rvalues$obs[rvalues$obs_cat]
    
    #group.by <- list(rvalues$h5ad[[1]]$obs[input$grouping_1][,,drop=TRUE])
    selectInput(
      inputId = 'de_anno_sel',
      label = 'Select the observations you want to compare!',
      choices = annotation_options,
      multiple = FALSE)
  })
  
  ## Baseline group for differential expression comparison
  output$de_group_1_list <- renderUI({
    req(input$assay_1,input$de_anno_sel)
    
    anno_choices <- unique(rvalues$h5ad[[1]]$obs[c(input$de_anno_sel)][,1])
    
    selectInput(
      inputId = 'de_group_1_sel',
      label = 'Select the annotation to use as baseline!',
      choices = anno_choices,
      multiple = FALSE)  
  })
  
  ## Group to compare against for differential expression comparison
  output$de_group_2_list <- renderUI({
    req(input$assay_1,input$de_anno_sel, input$de_group_1_sel)
    
    anno_choices_2 <- unique(rvalues$h5ad[[1]]$obs[c(input$de_anno_sel)][,1])
    anno_choices_2 <- setdiff(anno_choices_2,input$de_group_1_sel)
    
    selectInput(
      inputId = 'de_group_2_sel',
      label = 'Select the annotation to compare against!',
      choices = anno_choices_2,
      multiple = FALSE)  
  })
  
  ## Data frame containin expression for cells of DE group1
  de_analysis_group1 <- reactive({
    df_annos <- rvalues$h5ad[[1]]$obs[c(input$de_anno_sel)][,1]
    df_annos <- data.frame("anno" = df_annos)
    df_annos$cell_ids <- 1:nrow(df_annos)
    df_cells <- subset(df_annos,anno == input$de_group_1_sel)
    if(rvalues$raw_dtype == "NULL" | rvalues$raw_dtype == "counts"){ ## Check if raw exists
      df_cells_exp_anno <- rvalues$h5ad[[1]]$X[df_cells$cell_ids,]
    }else if(rvalues$raw_dtype == "normalized"){ ## Only if the other two conditions fail, use raw values to calculate differential expression
      df_cells_exp_anno <- rvalues$h5ad[[1]]$raw$X[df_cells$cell_ids,] 
    }
    df_cells_exp_anno
  })
  
  ## Data frame containing expression for cells of DE group2
  de_analysis_group2 <- reactive({
    df_annos <- rvalues$h5ad[[1]]$obs[c(input$de_anno_sel)][,1]
    df_annos <- data.frame("anno" = df_annos)
    df_annos$cell_ids <- 1:nrow(df_annos)
    df_cells <- subset(df_annos,anno == input$de_group_2_sel)
    if(rvalues$raw_dtype == "NULL" | rvalues$raw_dtype == "counts"){ ## Check if raw exists
      df_cells_exp_anno <- rvalues$h5ad[[1]]$X[df_cells$cell_ids,]
    }else if(rvalues$raw_dtype == "normalized"){ ## Only if the other two conditions fail, use raw values to calculate differential expression
      df_cells_exp_anno <- rvalues$h5ad[[1]]$raw$X[df_cells$cell_ids,] 
    }
    df_cells_exp_anno
  })
  
  ## Event when user clicks button to compare differential expression
  observeEvent(input$run_de_analysis,{
    
    de_res <- reactive({
      de_group1 <- isolate(de_analysis_group1())
      de_group2 <- isolate(de_analysis_group2())
      
      ## Merge both expression tables
      matrix_tr <- t(rbind(de_group1,de_group2))
      rownames(matrix_tr) <- rvalues$features
      
      ## Create a feature vector for both groups
      feature_vec_1 <- replicate(nrow(de_group1),isolate(input$de_group_1_sel))
      feature_vec_2 <- replicate(nrow(de_group2),isolate(input$de_group_2_sel))
      feature_vec <- c(feature_vec_1,feature_vec_2)
      
      de_res <- wilcoxauc(matrix_tr, feature_vec)
      de_res <- de_res %>%
        arrange(desc(auc)) %>%
        subset(group == isolate(input$de_group_2_sel)) %>%
        dplyr::select(-c(pct_in,pct_out,group,statistic))
      
      de_res
    })
    
    avg_exp <- reactive({
      req(de_res())
      de_group1 <- isolate(de_analysis_group1())
      de_group2 <- isolate(de_analysis_group2())
      avg_exp_group1 <- colMeans(de_group1)
      avg_exp_group2 <- colMeans(de_group2)
      avg_exp_df <- data.frame("group1" = avg_exp_group1,
                               "group2" = avg_exp_group2,
                               "gene" = rvalues$features)

      avg_exp_df <- avg_exp_df %>%
        mutate("significant" = if_else(gene %in% subset(de_res(),padj < 0.05)$feature,"yes","no"))
      avg_exp_df
    })
    
    output$de_res_table <- renderReactable({
      de_res_tbl_view <- de_res() 
      
      reactable(isolate(de_res_tbl_view),
                sortable = TRUE,
                searchable = TRUE,
                selection = "single")
    })
    
    selected_de <- reactive(getReactableState("de_res_table", "selected"))
    selected_de_gene <- reactive({
      isolate(de_res())[selected_de(),]$feature
    })
    
    ## Violin plot for differential expression
    output$de_volcano_plot <- renderPlot({
      req(de_res())
      de_res_tbl <- de_res() %>%
        mutate("padj" = if_else(padj == 0,1e-300,padj))
      
      volcano_plot <- ggplot(isolate(de_res_tbl),aes(logFC,-log10(padj))) +
        geom_point() +
        bbc_style() +
        theme(axis.title = element_text(size = 16, face = "bold")) +
        labs(title = "Volcano plot",
             x = "log2FC",
             y = "-log10(padj)") +
        scale_y_continuous(limits = c(0,max(-log10(de_res_tbl$padj))+ 30))
      if(!is.null(selected_de())){
        volcano_plot <- volcano_plot + 
          geom_point(data = subset(isolate(de_res_tbl), feature == selected_de_gene()),
                     size = 5, fill = "red",color=  "black",pch = 21)
      }
      volcano_plot
    })
    
    ## Correlation plot for differential expression
    output$de_avg_exp_plot <- renderPlot({
      group1_name <- isolate(input$de_group_1_sel)
      group2_name <- isolate(input$de_group_2_sel)
      avg_exp_plot <- ggplot(isolate(avg_exp()),aes(group1,group2)) +
        geom_point() +
        geom_abline(linetype = 2) +
        labs(title = "Average expression",
             x = group1_name,
             y = group2_name) +
        bbc_style() +
        theme(axis.title = element_text(size = 16, face = "bold"))
      
      if(!is.null(selected_de())){
        avg_exp_plot <- avg_exp_plot + 
          geom_point(data = subset(isolate(avg_exp()), gene == selected_de_gene()),
                     size = 5, fill = "red",color=  "black",pch = 21)
      }
      avg_exp_plot
    })
    
    de_violin_data <- reactive({

      ## Merge both expression tables
      de_group_1 <- isolate(de_analysis_group1())
      group1_name <- isolate(input$de_group_1_sel)
      de_group_2 <- isolate(de_analysis_group2())
      group2_name <- isolate(input$de_group_2_sel)
      
      matrix <- rbind(de_group_1,de_group_2)
      colnames(matrix) <- rvalues$features
      matrix_sub <- matrix[,selected_de_gene()]
      
      ## Create a feature vector for both groups
      feature_vec_1 <- replicate(nrow(de_group_1),group1_name)
      feature_vec_2 <- replicate(nrow(de_group_2),group2_name)
      feature_vec <- c(feature_vec_1,feature_vec_2)
      
      final_df <- data.frame("exp" = matrix_sub,
                             "group" = feature_vec)
      
      ## set order of the two groups, such that group1 is always first
      final_df$group <- factor(final_df$group,
                               levels = c(group1_name,group2_name))
        
      final_df
    })
    

    ## Correlation plot for differential expression
    output$de_violin_plot <- renderPlot({
      req(selected_de_gene())
      ggplot(isolate(de_violin_data()),aes(group,exp, fill = group)) +
        geom_violin() +
        stat_summary(fun=mean, geom="point", size=5, color = "black") +
        scale_fill_manual(values = c("#4682B4","#B47846")) +
        bbc_style() +
        theme(axis.title = element_text(size = 16, face = "bold")) +
        theme(legend.position = "none") +
        labs(x = "Feature",
             y = "Gene expression",
             title = selected_de_gene())
    })
    
    
  })
  
} # server end
