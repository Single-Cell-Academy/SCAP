#---------------------- Server ---------------------------#

server <- function(input, output, session){
  
  # set working dir and source required functions
  options(shiny.maxRequestSize=500*1024^2)
  setwd("/Users/jsjoyal/Desktop/SCAP/")
  volumes <- c(data = './test_data/')
  source("./R/SCAP_functions.R")
  
  # import data
  output$data_used <- renderText({
    if(!is.null(input$assay_1) & !is.na(input$directory[[1]][2])){
      return(paste0("Chosen data: ", input$directory[[1]][2]))
    }else{
      return(paste0("Please select your data..."))
    }
  })
  
  shinyDirChoose(input, "directory", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  # connect to .loom files in chosen dir
  data <- NULL
  data_update <- reactive({
    req(input$directory)
    path <- parseDirPath(volumes, input$directory)
    if(identical(path, character(0))) return(NULL)
    loom_files <- paste0(path,"/",list.files(path))
    assays <- sub(".loom","",sub(paste0(".*/"),"",loom_files))
    data <<- list()
    for(i in 1:length(assays)){
      data[[i]] <<- connect(loom_files[i], mode = "r+")
    }
    if(is.null(data)) return(NULL)
    names(data) <<- assays
    return(sample(1:10^9,1))
  })
  
  ###==============// MAIN TAB //==============####
  
  resolution.trigger <- makeReactiveTrigger()
  
  # global variables used by this tab
  grouping_options_old <- NULL
  
  observeEvent(input$low_res_1,{
    resolution.trigger$trigger()
  })
  
  #-- select input for Assay --#
  output$assay_1 <- renderUI({
    req(data_update())
    if(is.null(data)){
      return(NULL)
    }
    selectInput('assay_1', "Select Assay", choices = names(data), selected = ifelse(any(names(data)=="RNA"),yes = "RNA",no = names(data)[1]))
  })
  
  #-- select how to group cells --#
  output$grouping_1 <- renderUI({
    req(data_update())
    req(input$assay_1)
    assay <- input$assay_1
    cols <- names(data[[assay]]$col.attrs)
    grouping_options <- cols[grep('_meta_data$',cols)]
    if(any(grepl(paste0(tolower(assay),"_clusters"),sub("_meta_data","",grouping_options), fixed = TRUE))){
      sel <- paste0(tolower(assay),"_clusters")
    }else{
      sel <- sub("_meta_data","",grouping_options[1])
    }
    grouping_options_old <<- grouping_options
    
    selectInput(inputId = 'grouping_1', label = 'Group By', choices = sub("_meta_data$","",grouping_options), selected = sel, multiple = FALSE)
  })
  
  #-- select the clustering method --#
  output$reduction_1 <- renderUI({
    req(data_update())
    req(input$assay_1)
    assay <- input$assay_1
    options <- names(data[[assay]]$col.attrs)
    l.assay <- tolower(assay)
    reductions <- unique(sub("_._reduction$","",options[which(grepl('_reduction$',options)&grepl(l.assay,options))]))
    selectInput(inputId = 'reduction_1', 'Choose Clustering Method', choices = as.list(reductions), selected = ifelse(test = any(reductions==paste0('tsne_',l.assay,'_2d')), yes = paste0('tsne_',l.assay,'_2d'), no = reductions[1]))
  })
  
  #-- select the featutres for the feature plot --#
  output$featureplot_1_feature.select<- renderUI({
    req(data_update())
    req(input$assay_1)
    assay <- input$assay_1
    selectInput(inputId = 'featureplot_1_feature.select', label = 'Select a Feature to Visualize on the Feature Plot', choices = data[[assay]][['row_attrs/features']][], selected = data[[assay]][['row_attrs/features']][1], multiple = FALSE)
  })
  
  #-- select the featutres for the dot plot --#
  output$dotplot_1_feature_select<- renderUI({
    req(data_update())
    req(input$assay_1)
    assay <- input$assay_1
    selectInput(inputId = 'dotplot_1_feature_select', label = 'Choose Features to Visualize on Dot Plot', choices = data[[assay]][['row_attrs/features']][], selected = data[[assay]][['row_attrs/features']][3], multiple = TRUE)
  })
  
  #-- dimensional reduction plot coloured by cell groups --#
  output$dimplot_1 <- renderPlotly({
    req(data_update())
    req(input$grouping_1)
    req(input$reduction_1)
    resolution.trigger$depend()
    dimPlotlyOutput(assay.in = input$assay_1, reduc.in = input$reduction_1, group.by = input$grouping_1, annot_panel = "", low.res = 'no',data = data)
  })
  
  #-- dimensional reduction plot coloured by feature expression --#
  output$featureplot_1 <- renderPlotly({
    req(data_update())
    req(input$featureplot_1_feature.select)
    resolution.trigger$depend()
    featurePlotlyOutput(assay.in = input$assay_1, reduc.in = input$reduction_1, group.by = input$grouping_1, feature.in = input$featureplot_1_feature.select, low.res = 'no', data = data)
  })
  
  #-- dot plot for selected feature expression --#
  output$dotplot_1 <- renderPlot({
    req(data_update())
    req(input$dotplot_1_feature_select)
    
    group.by <- input$grouping_1
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
    p <- dotPlot(data = data[[assay]], assay = assay, features = input$dotplot_1_feature_select, group.by = group.by)
    if(is.null(p)) return(NULL)
    if(identical(class(p),"shiny.tag")) return(NULL)
    
    p <- p + guides(color = guide_colourbar(order = 1, title = "Average Expression"), size = guide_legend(order = 2, title = "Percent Expressed")) + theme_few()
    
    return(p)
  })
  
  ###==============// ANNOTATION TAB //==============####
  
  observeEvent(input$low_res_2,{
    resolution.trigger$trigger()
  })
  
  #-- global variables used by this tab --#
  tmp_annotations <- NULL
  meta_options_old <- NULL
  cells <- NULL
  
  #-- reactive triggers used by this tab --#
  # triggered when adding custom annotations
  tmp.trigger <- makeReactiveTrigger()
  # triggered when saving custom annotations
  annot.trigger <- makeReactiveTrigger()
  
  #-- select input for Assay --#
  output$assay_2 <- renderUI({
    req(data_update())
    selectInput('assay_2', "Select Assay", choices = names(data), selected = ifelse(any(names(data)=="RNA"),yes = "RNA",no = names(data)[1]))
  })
  
  #-- initialize/reset tmp_annotations --#
  observeEvent(input$assay_2,{
    annot.trigger$depend()
    req(data_update())
    tmp_annotations <<- rep("unlabled",data[[input$assay_2]]$shape[2])
    names(tmp_annotations) <<- data[[input$assay_2]]$col.attrs$CellID[]
  })
  
  #-- select how to group cells --#
  output$grouping_2 <- renderUI({
    annot.trigger$depend()
    req(data_update())
    req(input$assay_2)
    assay <- input$assay_2
    cols <- names(data[[assay]]$col.attrs)
    meta_options <- cols[grep('_meta_data$',cols)]
    if(!is.null(meta_options_old) & (length(meta_options)!=length(meta_options_old))){
      sel <- sub("_meta_data$", "", meta_options[!meta_options%in%meta_options_old])
    }else if(any(grepl(paste0(tolower(assay),"_clusters"),sub("_meta_data","",meta_options), fixed = TRUE))){
      sel <- paste0(tolower(assay),"_clusters")
    }else{
      sel <- sub("_meta_data","",meta_options[1])
    }
    meta_options_old <<- meta_options
    
    selectInput(inputId = 'grouping_2', label = 'Group By', choices = sub("_meta_data$","",meta_options), selected = sel, multiple = FALSE)
  })
  
  #-- select the clustering method --#
  output$reduction_2 <- renderUI({
    req(data_update())
    req(input$assay_2)
    assay <- input$assay_2
    options <- names(data[[assay]]$col.attrs)[grep("_2d",names(data[[assay]]$col.attrs))]
    l.assay <- tolower(assay)
    reductions <- unique(sub("_._reduction$","",options[which(grepl('_reduction$',options)&grepl(l.assay,options))]))
    selectInput(inputId = 'reduction_2', 'Choose Clustering Method', choices = as.list(reductions), selected = ifelse(test = any(reductions==paste0('tsne_',l.assay,'_2d')), yes = paste0('tsne_',l.assay,'_2d'), no = reductions[1]))
  })
  
  #-- select the featutres for the feature plot --#
  output$featureplot_2_feature_select<- renderUI({
    req(data_update())
    req(input$assay_2)
    assay <- input$assay_2
    selectInput(inputId = 'featureplot_2_feature_select', label = 'Select a Feature to Visualize on the Feature Plot', choices = data[[assay]][['row_attrs/features']][], selected = data[[assay]][['row_attrs/features']][1], multiple = FALSE)
  })
  
  #-- display cells that were selected in a table --#
  output$selected_cells <- renderTable({
    req(data_update())
    req(input$assay_2)
    selected_cells <- as.data.frame(event_data("plotly_selected")$key, stringsAsFactors = FALSE)
    if (is.null(selected_cells) | ncol(selected_cells) == 0 | nrow(selected_cells) == 0){
      cells <<- NULL
      return("Click and drag on either plot to select cells for marker identification (i.e., select/lasso) (double-click to clear selection)")
    }else if(ncol(selected_cells)>1){
      selected_cells <- as.data.frame(as.character(selected_cells[1,]),stringsAsFactors = F)
    }
    colnames(selected_cells) <- ""
    cells <<- selected_cells
    selected_cells
  })
  
  #-- dimensional reduction plot coloured by cell groups --#
  output$dimplot_2 <- renderPlotly({
    resolution.trigger$depend()
    req(data_update())
    req(input$grouping_2)
    tmp.trigger$depend()
    annot.trigger$depend()
    dimPlotlyOutput(assay.in = input$assay_2, reduc.in = input$reduction_2, group.by = input$grouping_2, annot_panel = input$annot_panel, tmp_annotations = tmp_annotations, low.res = input$low_res_2, data = data)
  })
  
  #-- dimensional reduction plot coloured by feature expression --#
  output$featureplot_2 <- renderPlotly({
    resolution.trigger$depend()
    req(data_update())
    req(input$featureplot_2_feature_select)
    tmp.trigger$depend()
    annot.trigger$depend()
    featurePlotlyOutput(assay.in = input$assay_2, reduc.in = input$reduction_2, group.by = input$grouping_2, feature.in = input$featureplot_2_feature_select, low.res = input$low_res_2, data = data)
  })
  
  #-- Find the markers for the 1) selected cells compared to all other cells or
  #-- 2) the markers that define each group. Depending on if cells are selected or not
  observeEvent(input$find.markers, ignoreNULL = FALSE, ignoreInit = FALSE, handlerExpr = {
    output$markers <- renderTable({
      req(data_update())
      req(input$grouping_2)
      grouping <- paste0(input$grouping_2,"_meta_data")
      if(!any(names(data[[input$assay_2]]$col.attrs) == grouping)){
        return(NULL)
      }
      
      m <- t(data[[input$assay_2]]$matrix[,])
      rownames(m) <- data[[input$assay_2]]$row.attrs$features[]
      colnames(m) <- data[[input$assay_2]]$col.attrs$CellID[]
      
      y <- data[[input$assay_2]][[paste0('col_attrs/',grouping)]][]
      cols <- unlist(lapply(cells[,1], function(x){
        grep(x,data[[input$assay_2]]$col.attrs$CellID[])
      }))
      
      if(!is.null(cells))
        y[cols] <- 'Selected'
      
      t <- as.data.frame(top_markers(wilcoxauc(X = m, y = y), n = 10))[1:10,]
      if(!is.null(cells)){
        t <- as.data.frame(t[,which(colnames(t)=="Selected")],stringsAsFactors=FALSE)
        names(t) <- "Markers for Selected Cells"
      }else{
        t <- t[,mixedsort(colnames(t))]
        index <- 1:ncol(t)
        rank <- grep('rank',colnames(t))
        index <- index[-rank]
        t <- t[,c(rank, index)]
      }
      t
    }, hover = TRUE, width = '100%',align = 'c')
  })
  
  output$annotation_title <- renderText({
    req(data_update())
    annot.trigger$depend()
    req(input$grouping_2)
    return(paste0(input$grouping_2))
  })
  
  #-- display text boxes for user to enter new IDs for the meta data --#
  output$annotations <- renderUI({
    annot.trigger$depend()
    req(data_update())
    req(input$assay_2)
    req(input$grouping_2)
    group.by <- paste0(input$grouping_2,'_meta_data')
    if(!group.by%in%names(data[[input$assay_2]]$col.attrs)){
      return(NULL)
    }
    names <- unique(data[[input$assay_2]][[paste0('col_attrs/',group.by)]][])
    max_names <- 50
    if(length(names)>max_names){
      showNotification(paste0('Warning: Adding annotations to a grouping of this size (',length(names),') is not supported.'), type = 'warning')
      return(NULL)
    }
    names <- mixedsort(names)
    id <- "my_annot"
    lapply(1:length(names), function(x){
      textInput(inputId = paste0(names[x],"_",id), label = names[x], placeholder = names[x])
    })
  })
  
  #-- add user defined annotations to selected cells --#
  observeEvent(input$add_to_tmp, {
    if(is.null(cells)){
      showNotification("Please select cells from either plot first", type = 'warning')
    }else if(is.null(input$custom_name) | input$custom_name == ""){
      showNotification("Please enter an annotation name", type = 'warning')
    }else{
      tmp_annotations[which(names(tmp_annotations)%in%cells[,1])] <<- input$custom_name
    }
    tmp.trigger$trigger()
  })

  #-- save user defined annotations --#
  observeEvent(input$set_cell_types, ignoreNULL = FALSE, ignoreInit = FALSE, handlerExpr = {
    req(data_update())
    req(input$grouping_2)
    group.by <- paste0(input$grouping_2,'_meta_data')
    id <- "my_annot"
    if(is.null(input$new_scheme_name) || input$new_scheme_name == "" || input$new_scheme_name == " "){
      showNotification("Please name the annotation scheme", type = 'warning')
    }else if(input$annot_panel == 'cell_annotation_cluster'){
      names <- data[[input$assay_2]][[paste0('col_attrs/',group.by)]][]
      new <- names
      names <- unique(names)
      for(i in 1:length(names)){
        new[which(new == names[i])] <- input[[paste0(names[i],"_",id)]]
        if(is.null(input[[paste0(names[i],"_",id)]]) | input[[paste0(names[i],"_",id)]] == ""){
          showNotification("Names must be provided for each group", type = 'warning')
          return()
        }
      }
      new <- list(new)
      names(new) <- input$new_scheme_name
      if(!grepl("_meta_data$",names(new))) names(new) <- paste0(names(new),"_meta_data")
      
      data[[input$assay_2]]$add.col.attribute(attributes = new, overwrite = TRUE)
      annot.trigger$trigger()
    }else{
      if(all(tmp_annotations=='unlabled')){
        showNotification("No annotations have been added", type = 'warning')
      }else if(is.null(input$new_scheme_name) | input$new_scheme_name==""){
        showNotification("Please enter an name for the annotation scheme", type = 'warning')
      }else{
        new <- list(tmp_annotations)
        names(new) <- input$new_scheme_name
        if(!grepl("_meta_data$",names(new))){
          names(new) <- paste0(names(new),"_meta_data")
        }
        data[[input$assay_2]]$add.col.attribute(attributes = new, overwrite = TRUE)
        annot.trigger$trigger()
      }
    }
  })
  
  ###==============// CUSTOM META DATA TAB //==============####
  
  meta.trigger <- makeReactiveTrigger()

  order <- NULL
  
  output$assay_3 <- renderUI({
    req(data_update())
    selectInput('assay_3', "Select Assay", choices = names(data), selected = ifelse(any(names(data)=="RNA"),yes = "RNA",no = names(data)[1]))
  })
  
  output$meta_group_checkbox <- renderUI({
    req(data_update())
    req(input$assay_3)
    choices <- names(data[[input$assay_3]]$col.attrs)[which(grepl("_meta_data$", names(data[[input$assay_3]]$col.attrs)))]
    choices <- sub('_meta_data','',choices)
    checkboxGroupInput(inputId = 'meta_group_checkbox', label = 'Meta Data', choices = choices)
  })
  
  output$example_meta_group_name <- renderText({
    req(data_update())
    meta.trigger$depend()
    if(is.null(order) || is.na(order) || length(order) == 0){
      return(NULL)
    }
    paste(sub('_meta_data$','',order),collapse=" + ")
  })
  
  output$example_meta_group <- renderTable({
    req(data_update())
    req(input$assay_3)
    meta.trigger$depend()
    if(is.null(order) || is.na(order) || length(order) == 0){
      return(NULL)
    }
    ex <- data[[input$assay_3]]$get.attribute.df(attributes = order)
    ex_df <- as.data.frame(apply(X = ex, MARGIN = 1, FUN = paste, collapse = '_'), stringsAsFactors = FALSE)[sample(x = 1:nrow(ex),size = 20,replace = FALSE),,drop=F]
    colnames(ex_df) <- 'Examples of chosen meta data combination'
    return(ex_df)
  }, hover = TRUE, width = '100%',align = 'c')
  
  observeEvent(input$meta_group_checkbox, ignoreNULL = F,{
    req(data_update())
    if(is.null(new)){
      order <<- NULL
    }else if(!is.null(input$meta_group_checkbox)){
      new <- paste0(input$meta_group_checkbox,'_meta_data')
      if(length(new)>length(order) | is.null(order)){
        new <- new[!new%in%order]
        order <<- c(order, new)
      }else if(identical(new, character(0))){
        order <<- NULL
      }else{
        order <<- order[-which(!order%in%new)]
      }
    }else{
      order <<- NULL
    }
    meta.trigger$trigger()
  })
}
