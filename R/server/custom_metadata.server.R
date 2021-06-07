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