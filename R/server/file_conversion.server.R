###==============// FILE CONVERSION TAB //==============####

shinyFileChoose(input, "file_1", roots = volumes, session = session)
shinyDirChoose(input, "out_dir_2", roots = volumes, session = session, restrictions = system.file(package = "base"))

output$chosen_in_file <- renderText(parseFilePaths(selection = input$file_1, roots = volumes)$datapath)
output$chosen_out_directory <- renderText(parseDirPath(roots = volumes, selection = input$out_dir_2))
output$chosen_out_file <- renderText(input$out_file_2)

observeEvent(input$scap_convert, ignoreInit = T, ignoreNULL = T, {
  file_1 <- parseFilePaths(selection = input$file_1, roots = volumes)$datapath
  if(identical(file_1,character(0))){
    return(NULL)
  }
  if(grepl('\\.rds|\\.h5ad|\\.loom', file_1, ignore.case = TRUE) == FALSE){
    message('Error: Unsupported File Format')
    return(NULL)
  }
  file_2 <- paste0(parseDirPath(roots = volumes, selection = input$out_dir_2), "/", input$out_file_2)
  showModal(modalDialog(p("Converting, Please Wait..."), title = "This window will close after conversion is complete"), session = getDefaultReactiveDomain())
  SCAP_Convert(input$from, input$to, file_1, file_2)
  removeModal(session = getDefaultReactiveDomain())
})