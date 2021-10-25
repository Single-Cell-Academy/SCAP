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
                  reduc.in = rvalues$reductions[[input$reduction_2]][,c(1,2)], 
                  group.by = group.by, 
                  annot_panel = input$annot_panel, 
                  tmp_annotations = rvalues$tmp_annotations, 
                  low.res = 'yes',
                  hide.legend = input$hidelegend_2)
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
                      reduc.in = rvalues$reductions[[input$reduction_2]][,c(1,2)], 
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
      
      scanpy$tl$rank_genes_groups(rvalues$h5ad[[1]], groupby = 'scap_find_markers_groups', use_raw = use_raw(rvalues$h5ad[[1]]), method = 'wilcoxon') # find DEGs
      
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
      names <- as.character(rvalues$h5ad[[1]]$obs[group.by][,,drop = TRUE])
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