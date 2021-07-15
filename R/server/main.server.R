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
    sel <- rvalues$obs[grep("seurat_clusters", rvalues$obs, ignore.case = TRUE)]
  }else if(any(grepl(paste0(tolower(assay),"_clusters"), rvalues$obs, ignore.case = TRUE))){
    sel <- rvalues$obs[grep(paste0(tolower(assay),"_clusters"), rvalues$obs, ignore.case = TRUE)]
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
    if(length(unique(rvalues$h5ad[[1]]$obs[x][,1,drop=TRUE]))>1){
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
  group.by <- list(rvalues$h5ad[[1]]$obs[input$grouping_1][,1,drop=TRUE])
  names(group.by) <- input$grouping_1
  names(group.by[[1]]) <- rvalues$cell_ids
  cat <- rvalues$obs_cat[which(rvalues$obs == input$grouping_1)]
  if(cat){
    dimPlotlyOutput(assay.in = input$assay_1, 
                    reduc.in = rvalues$reductions[[input$reduction_1]][,c(1,2)], 
                    group.by = group.by, 
                    annot_panel = "", 
                    low.res = 'yes')
  }else{
    feature.in <- rvalues$h5ad[[1]]$obs[input$grouping_1][,1,drop=FALSE]
    colnames(feature.in) <- input$grouping_1
    rownames(feature.in) <- rownames(rvalues$reductions[[input$reduction_1]])
    featurePlotlyOutput(assay.in = input$assay_1,
                        reduc.in = rvalues$reductions[[input$reduction_1]][,c(1,2)],
                        group.by = group.by,
                        feature.in = feature.in,
                        low.res = 'yes')
  }
})

#-- dimensional reduction plot coloured by feature expression --#
output$featureplot_1 <- renderPlotly({
  req(input$grouping_1, input$reduction_1, input$featureplot_1_feature_select)
  group.by <- list(rvalues$h5ad[[1]]$obs[input$grouping_1][,1,drop=TRUE])
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
                        reduc.in = rvalues$reductions[[input$reduction_1]][,c(1,2)],
                        group.by = group.by,
                        feature.in = feature.in,
                        low.res = 'yes')
  }else{
    featurePlotlyOutput_nebulosa(assay.in = input$assay_1,
                                 reduc.in = rvalues$reductions[[input$reduction_1]][,c(1,2)],
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
  data.features$id <- rvalues$h5ad[[1]]$obs[input$grouping_1][,1,drop=TRUE]
  
  if(input$do_split == 'yes' & !is.null(input$split_by)){
    splits = list(rvalues$h5ad[[1]]$obs[input$split_by][,1,drop=TRUE])
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