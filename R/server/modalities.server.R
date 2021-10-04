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
  if((length(data[[1]]$obs_names) != length(rvalues$h5ad[[1]]$obs_names)) || (all(as.character(data[[1]]$obs_names) == as.character(rvalues$h5ad[[1]]$obs_names)) == FALSE)){
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
  
  ## Determine what data is likely stored in .raw
  if(is.null(data[[1]]$raw)){ ## Check if raw exists
    rvalues_mod$raw_dtype <- "NULL"
  }else if(sum(rvalues$h5ad[[1]]$raw$X[1,]) %% 1 == 0){ ## Check whether raw contains un-normalized data or normalized data
    rvalues_mod$raw_dtype <- "counts"
  }else{ ## Only if the other two conditions fail, use raw values to calculate differential expression
    rvalues_mod$raw_dtype <- "normalized"
  }
  
})

#-- select input for Assay --#
output$assay_mod <- renderUI({
  if(is.null(rvalues_mod$h5ad)){
    return(NULL)
  }else{
    return(selectInput('assay_mod', "Select Assay", choices = names(rvalues_mod$h5ad), selected = ifelse(any(names(rvalues_mod$h5ad)=="RNA"), yes = "RNA", no = names(rvalues_mod$h5ad)[1])))
  }
})

#-- Display name of data --#
output$data_used_mod <- renderUI({
  if(is.null(rvalues_mod$h5ad)){
    h3('Select data')
  }else{
    path <- sub(".*\\/", "", parseFilePaths(selection = input$h5ad_in_mod, roots = volumes)$datapath)
    h3(paste0("Chosen data: ", path))
  }
})

output$grouping_mod <- renderUI({
  req(rvalues_mod$h5ad, input$assay_mod, rvalues_mod$obs)
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
  req(rvalues_mod$h5ad, input$assay_mod)
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
  req(rvalues_mod$h5ad, input$assay_mod)
  assay <- input$assay_mod
  selectInput(inputId = 'featureplot_mod_feature_select', 
              label = 'Select a Feature to Visualize on the Feature Plot', 
              choices = rvalues_mod$features, 
              selected = rvalues_mod$features[1], 
              multiple = ifelse(input$nebulosa_on == "yes", TRUE, FALSE))
})

output$ridgeplot_mod_feature_select <- renderUI({
  req(rvalues_mod$h5ad, input$assay_mod)
  assay <- input$assay_mod
  selectInput(inputId = 'ridgeplot_mod_feature_select', 
              label = 'Select Features to Visualize on the Ridge Plot', 
              choices = rvalues_mod$features, 
              selected = rvalues_mod$features[1], 
              multiple = TRUE)
})

#-- dimensional reduction plot coloured by cell groups --#
output$dimplot_mod <- renderPlotly({
  req(rvalues_mod$h5ad, input$grouping_mod, input$reduction_mod)
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
                        reduc.in = rvalues_mod$reductions[[input$reduction_mod]][,c(1,2)],
                        group.by = group.by,
                        feature.in = feature.in,
                        low.res = 'yes')
  }
})

#-- dimensional reduction plot coloured by feature expression --#
output$featureplot_mod <- renderPlotly({
  req(rvalues_mod$h5ad, input$grouping_mod, input$reduction_mod, input$featureplot_mod_feature_select)
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
                        reduc.in = rvalues_mod$reductions[[input$reduction_mod]][,c(1,2)],
                        group.by = group.by,
                        feature.in = feature.in,
                        low.res = 'yes')
  }else{
    featurePlotlyOutput_nebulosa(assay.in = input$assay_mod,
                                 reduc.in = rvalues_mod$reductions[[input$reduction_mod]][,c(1,2)],
                                 group.by = as.data.frame(group.by),
                                 feature.in = feature.in,
                                 low.res = 'yes')
  }
})

## output variable to hold type of user selected grouping for main panel
output$grouping_mod_type <- reactive({
  req(rvalues_mod$h5ad, input$grouping_mod)
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
  req(rvalues_mod$h5ad, input$assay_mod)
  selectInput(inputId = 'crispr_feature_1_sel', 
              label = 'Select feature #1 you want to compare!', 
              choices = rvalues_mod$features, 
              selected = rvalues_mod$features[1], 
              multiple = FALSE)
})

output$crispr_feature_1_slider <- renderUI({ ## Feature 1 selected for CRISPR
  req(rvalues_mod$h5ad, input$assay_mod)
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
  req(rvalues_mod$h5ad, input$assay_mod)
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
  req(rvalues_mod$h5ad, input$assay_mod)
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
  req(rvalues_mod$h5ad, crispr_feature_1_cells())
  cells_pass <- subset(crispr_feature_1_cells(),pass_exp_thr == "yes")
  res_string <- paste("There are ",nrow(cells_pass)," cells above your selected threshold for feature:",
                      input$crispr_feature_1_sel,sep="")
  res_string
})


#### CRISPR feature 2 parameters
output$crispr_feature_2 <- renderUI({ ## Feature 1 selected for CRISPR
  req(rvalues_mod$h5ad, input$assay_mod)
  req(input$crispr_feature_1_sel)
  
  feature_options <- setdiff(rvalues_mod$features,input$crispr_feature_1_sel)
  
  selectInput(inputId = 'crispr_feature_2_sel', 
              label = 'Select feature #2 you want to compare!', 
              choices = feature_options, 
              selected = feature_options[1], 
              multiple = FALSE)
})

output$crispr_feature_2_slider <- renderUI({ ## Feature 1 selected for CRISPR
  req(rvalues_mod$h5ad, input$assay_mod)
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
  req(rvalues_mod$h5ad, input$assay_mod)
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
  req(rvalues_mod$h5ad, input$assay_mod)
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
  req(rvalues_mod$h5ad, crispr_feature_2_cells())
  cells_pass <- subset(crispr_feature_2_cells(),pass_exp_thr == "yes")
  res_string <- paste("There are ",nrow(cells_pass)," cells above your selected threshold for feature:",
                      input$crispr_feature_2_sel,sep="")
  res_string
})

#### Action button for calculating CRISPR DE
observeEvent(input$crispr_de_analysis,{
  req(rvalues_mod$h5ad)
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
  
  ## Calculate differential expression for RNA of both features
  mod_de_res <- reactive({
    de_group1 <- isolate(crispr_feature_1_cells_rna)
    de_group2 <- isolate(crispr_feature_2_cells_rna)
    
    ## Merge both expression tables
    matrix_tr <- t(rbind(de_group1,de_group2))
    rownames(matrix_tr) <- rownames(rvalues$h5ad[[1]]$var)
    
    ## Create a feature vector for both groups
    feature_vec_1 <- replicate(nrow(de_group1),isolate(input$crispr_feature_1_sel))
    feature_vec_2 <- replicate(nrow(de_group2),isolate(input$crispr_feature_2_sel))
    feature_vec <- c(feature_vec_1,feature_vec_2)
    
    de_res <- wilcoxauc(matrix_tr, feature_vec)
    de_res <- de_res %>%
      arrange(desc(logFC)) %>%
      subset(group == isolate(input$crispr_feature_2_sel)) %>%
      dplyr::select(-c(pct_in,pct_out,group,statistic))
    
    de_res
  })
  
  ## Table containing average RNA expression values for both features
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
    req(rvalues_mod$h5ad, merged_crispr_feature_avg_exp())
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
  
  selected <- reactive(getReactableState("crispr_de", "selected"))
  output$crispr_de <- renderReactable({
    req(mod_de_res())
    
    mod_de_res_tbl_view <- mod_de_res() 
    
    reactable(isolate(mod_de_res_tbl_view),
              sortable = TRUE,
              searchable = TRUE,
              selection = "single")
  })
  
  output$crispr_gene_vlnplot <- renderPlot({
    req(rvalues_mod$h5ad, mod_de_res())
    req(selected())
    
    gene_selected <- mod_de_res()[selected(),]$feature
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


## Correlation plot functions
output$corr_grouping <- renderUI({
  req(rvalues_mod$obs,input$assay_mod)
  assay <- input$assay_mod
  options <- rvalues_mod$obs
  if(any(grepl("seurat_clusters", options, ignore.case = FALSE))){
    sel <- rvalues$obs[grep("seurat_clusters", options)]
  }else if(any(grepl(paste0(tolower(assay),"_clusters"), options, ignore.case = TRUE))){
    sel <- paste0(tolower(assay),"_clusters")
  }else{
    sel <- options[1]
  }
  selectInput(inputId = 'corr_grouping', label = 'Group By', choices = options, selected = sel, multiple = FALSE)
})

output$corr_sub_grouping <- renderUI({
  req(input$corr_grouping)
  if(rvalues$obs_cat[which(rvalues$obs == input$corr_grouping)]){
    options <- c('All Options', levels(reorder_levels(rvalues$h5ad[[1]]$obs[[input$corr_grouping]])))
    selectInput(inputId = 'corr_sub_grouping', label = 'Show', choices = options, selected = options[1], multiple = FALSE)
  }else{
    sliderInput(inputId = 'corr_sub_grouping', label = 'Show', min = min(rvalues$h5ad[[1]]$obs[[input$corr_grouping]]), max = max(rvalues$h5ad[[1]]$obs[[input$corr_grouping]]), value = c(min(rvalues$h5ad[[1]]$obs[[input$corr_grouping]]), max(rvalues$h5ad[[1]]$obs[[input$corr_grouping]])))
  }
})

corr_df <- reactive({
  req(input$corr_fs_1, input$corr_fs_2, input$corr_grouping)
  
  df <- data.frame(x = rvalues$h5ad[[1]]$X[,which(rvalues$features == input$corr_fs_1), drop = TRUE], 
                   y = rvalues_mod$h5ad[[1]]$X[,which(rvalues_mod$features == input$corr_fs_2), drop = TRUE],
                   group = rvalues$h5ad[[1]]$obs[[input$corr_grouping]],
                   stringsAsFactors = FALSE
  )
  
  # add color pal
  if(rvalues$obs_cat[which(rvalues$obs == input$corr_grouping)]){
    df$col <- 'black'
    df$group <- reorder_levels(df$group)
    for(i in 1:length(levels(df$group))){
      df$col[which(df$group == levels(df$group)[i])] <- COLORPAL_DISCRETE[i]
    }
  }else{
    df <- df[order(df$group),]
    df$col <- COLORPAL_CONTINUOUS(nrow(df))
  }
  if(rvalues$obs_cat[which(rvalues$obs == input$corr_grouping)]){
    cat <- 1
  }else{
    cat <- 0
  }
  
  return(list(df = df, cat = cat))
})

corr_sub_df <- reactive({
  req(input$corr_sub_grouping)
  if(corr_df()$cat){
    if(identical(input$corr_sub_grouping, "All Options") == FALSE){
      df <- corr_df()$df[which(corr_df()$df$group %in% input$corr_sub_grouping),]
    }else{
      df <- corr_df()$df
    }
  }else{
    df <- corr_df()$df[which(corr_df()$df$group >= input$corr_sub_grouping[1] & corr_df()$df$group <= input$corr_sub_grouping[2]),]
  }
  
  flag <- (length(unique(df$x)) == 1 | length(unique(df$y)) == 1)
  if(flag){
    fit <- NULL
  }else{
    fit <- summary(lm(y~x, df))
  }
  return(list(df = df, fit = fit))
})

output$corr_fs_1 <- renderUI({
  req(input$assay_1)
  assay <- input$assay_1
  selectInput(inputId = 'corr_fs_1', 
              label = paste0('Select Feature 1 from ', assay), 
              choices = rvalues$features, 
              selected = rvalues$features[1], 
              multiple = FALSE)
})

output$corr_fs_2 <- renderUI({
  req(input$assay_mod)
  assay <- input$assay_mod
  selectInput(inputId = 'corr_fs_2', 
              label = paste0('Select Feature 2 from ', assay), 
              choices = rvalues_mod$features, 
              selected = rvalues_mod$features[1], 
              multiple = FALSE)
})

output$corr_plot <- renderPlot({
  req(corr_sub_df())
  p <- ggplot(data = corr_sub_df()$df, mapping = aes(x = x, y = y)) + 
    geom_point(color = corr_sub_df()$df$col) + 
    xlab(input$corr_fs_1) + 
    ylab(input$corr_fs_2) +
    theme_base() + 
    theme(plot.background = element_blank())
  if(is.null(corr_sub_df()$fit) == FALSE){
    p <- p + geom_abline(slope = corr_sub_df()$fit$coefficients[2,1], intercept = corr_sub_df()$fit$coefficients[1,1], color = "red")
  }
  return(p)
})

output$corr_stats <- renderUI({
  req(corr_sub_df())
  if(is.null(corr_sub_df()$fit)){
    HTML("One or more of the selected features has constant expression.")
  }else{
    HTML(paste(c(paste('slope:', round(corr_sub_df()$fit$coefficients[2,1],4)), paste('intercept:', round(corr_sub_df()$fit$coefficients[1,1],4)), paste('R2:', round(corr_sub_df()$fit$r.squared,4))), collapse = "<br/>"))
  }
})