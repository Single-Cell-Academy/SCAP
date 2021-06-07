## Scibet server elements

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