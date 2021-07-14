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