#' Shiny app server function
#'
#' @import cowplot
#' @import dplyr
#' @import ggplot2
#' @import ggthemes
#' @import gtools
#' @import hdf5r
#' @import loomR
#' @import Matrix
#' @import MODIS
#' @import plotly
#' @import presto
#' @import Seurat
#' @import shiny
#' @import shinycssloaders
#' @import shinyFiles
#' @import shinyjqui
#' @import shinythemes
#' @import rjson

library("cowplot")
library("dplyr")
library("ggplot2")
library("ggthemes")
library("gtools")
library("loomR")
library("hdf5r")
library("Matrix")
library("MODIS")
library("plotly")
library("presto")
library("Seurat")
library("rjson")
library("shiny")
library("shinycssloaders")
library("shinyFiles")
library("shinyjqui")
library("shinythemes")



server <- function(input, output, session){
  
  session$onSessionEnded(stopApp)
  
  # set working dir and source required functions
  options(shiny.maxRequestSize=500*1024^2)
  # currently requires manual changing to conform with users file system...should be fixed.
  # setwd("/Users/jsjoyal/Desktop/SCAP/")
  
  volumes <- c(projects = "/Users/jsjoyal/") #'/home/joel/SCAP/test_data/projects/'   #getVolumes()
  source("SCAP_functions.R")
  
  # import data
  output$data_used <- renderText({
    if(!is.null(input$assay_1) & !is.na(input$directory[[1]][2])){
      return(paste0("Chosen data: ", input$directory[[1]][2]))
    }else{
      return(paste0("Please select your data..."))
    }
  })
  
  shinyDirChoose(input, "directory", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  
  update.trigger <- makeReactiveTrigger()
  
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
    update.trigger$depend()
    return(sample(1:10^9,1))
  })
  
  
  # global variables
  m <- NULL # annotations tab: stores the normalized count matrix as a sparse matrix for speed
  tmp_annotations <- NULL # annotations tab: temporarily stores the custom annonations supplied by the user
  meta_options_old <- NULL # depricated
  cells <- NULL # annotations tab: stores a vector of cell IDs selected by the user (via plotly lasso)
  marker_tables <- list() # annotation tab: stores the tables of top 10 marker features for speed
  order <- NULL # custom metadata tab: stores the order of slected meta data by which order it was selected, rather than alpha-numerically 
  
  ###==============// MAIN TAB //==============####
  
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
    annot.trigger$depend()
    req(data_update())
    req(input$assay_1)
    assay <- input$assay_1
    cols <- names(data[[assay]]$col.attrs)
    grouping_options <- cols[grep('_meta_data$',cols)]
    df <- data[[assay]]$get.attribute.df(attributes = grouping_options)
    pass <- unlist(lapply(grouping_options, function(x){length(unique(df[,x]))<50}))
    grouping_options <- grouping_options[which(pass)]
    if(any(grepl(paste0(tolower(assay),"_clusters"),sub("_meta_data","",grouping_options), fixed = TRUE))){
      sel <- paste0(tolower(assay),"_clusters")
    }else{
      sel <- sub("_meta_data","",grouping_options[1])
    }
    selectInput(inputId = 'grouping_1', label = 'Group By', choices = sub("_meta_data$","",grouping_options), selected = sel, multiple = FALSE)
  })
  
  #-- select the clustering method --#
  output$reduction_1 <- renderUI({
    req(data_update())
    req(input$assay_1)
    assay <- input$assay_1
    options <- names(data[[assay]]$col.attrs)
    l.assay <- tolower(assay)
    reductions <- unique(sub("_._reduction$","",options[which(grepl('_reduction$',options))]))
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
  
  #-- visualize dotplot as a split dotplot --#
  output$do_split <- renderUI({
    req(data_update())
    req(input$assay_1)
    radioButtons(inputId = 'do_split', label = 'Split dot plot', choices = c('yes', 'no'), selected = 'no', inline = TRUE)
  })
  
  #-- select how to split the dot plot --#
  output$split_by <- renderUI({
    req(data_update())
    req(input$assay_1)
    choices <- unlist(lapply(names(data[[input$assay_1]]$col.attrs), function(x){
      if(length(unique(data[[input$assay_1]]$col.attrs[[x]][]))==2){
        return(x)
      }
    }))
    if(is.null(choices)){
      return("No conditions to split by.")
    }else{
      choices <- sub("_meta_data$","",choices[grepl('_meta_data$',choices)])
      return(radioButtons(inputId = 'split_by', label = 'Choose how to split the data', choices = choices, inline = FALSE))
    }
  })
  
  #-- dimensional reduction plot coloured by cell groups --#
  output$dimplot_1 <- renderPlotly({
    req(data_update())
    req(input$grouping_1)
    req(input$reduction_1)
    dimPlotlyOutput(assay.in = input$assay_1, reduc.in = input$reduction_1, group.by = input$grouping_1, annot_panel = "", low.res = 'yes',data = data)
  })
  
  #-- dimensional reduction plot coloured by feature expression --#
  output$featureplot_1 <- renderPlotly({
    req(data_update())
    req(input$featureplot_1_feature.select)
    featurePlotlyOutput(assay.in = input$assay_1, reduc.in = input$reduction_1, group.by = input$grouping_1, feature.in = input$featureplot_1_feature.select, low.res = 'yes', data = data)
  })
  
  #-- dot plot for selected feature expression --#
  dot_plot <- reactive({
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
    
    if(input$do_split == 'yes' & !is.null(input$split_by)){
      p <- split_dot_plot(data = data[[assay]], features = input$dotplot_1_feature_select, group.by = group.by, split.by = input$split_by)
      if(is.null(p)) return(NULL)
      if(identical(class(p),"shiny.tag")) return(NULL)
      return(p)
    }else if(input$do_split == 'no'){
      p <- dotPlot(data = data[[assay]], assay = assay, features = input$dotplot_1_feature_select, group.by = group.by)
      if(is.null(p)) return(NULL)
      if(identical(class(p),"shiny.tag")) return(NULL)
      
      p <- p + guides(color = guide_colourbar(order = 1, title = "Average Expression"), size = guide_legend(order = 2, title = "Percent Expressed")) + theme_few()
      
      return(p)
    }else{
      return(-1)
    }
  })
  
  output$dotplot_1 <- renderUI({
    req(data_update())
    req(dot_plot())
    req(input$dotplot_1_feature_select)
  
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
  
  #-- store tables of top 10 marker features --#
  #-- and store normalized count matrix in memory for speed --#
  observeEvent(data_update(), ignoreInit = TRUE, {
    showModal(modalDialog(p("Initializing data for better performance...This may take a little while."), title = paste0("Getting things set up for ", input$directory[[1]][2])), session = getDefaultReactiveDomain())
    marker_tables <<- list()
    assay_matrices <<- list()
    for(x in names(data)){
      assay_matrices[[x]] <<- t(Matrix(data = data[[x]]$matrix[,], sparse = TRUE))
      rownames(assay_matrices[[x]]) <<- data[[x]]$row.attrs$features[]
      colnames(assay_matrices[[x]]) <<- data[[x]]$col.attrs$CellID[]
    }
    cells <<- NULL
    removeModal(session = getDefaultReactiveDomain())
  })
  
  #-- reactive triggers used by this tab --#
  # triggered when adding custom annotations
  tmp.trigger <- makeReactiveTrigger()
  # triggered when saving custom annotations
  annot.trigger <- makeReactiveTrigger()
  #
  matrix.trigger <- makeReactiveTrigger()
  
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
    df <- data[[assay]]$get.attribute.df(attributes = meta_options)
    pass <- unlist(lapply(meta_options, function(x){length(unique(df[,x]))<50}))
    meta_options <- meta_options[which(pass)]
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
    reductions <- unique(sub("_._reduction$","",options[which(grepl('_reduction$',options))]))
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
      return("Click and drag on either plot to select cells for marker identification (i.e., select/lasso) (double-click on either plot to clear selection)")
    }else if(ncol(selected_cells)>1){
      selected_cells <- as.data.frame(as.character(selected_cells[1,]),stringsAsFactors = F)
    }
    colnames(selected_cells) <- ""
    cells <<- selected_cells
    selected_cells
  })
  
  #-- dimensional reduction plot coloured by cell groups --#
  output$dimplot_2 <- renderPlotly({
    req(data_update())
    req(input$grouping_2)
    tmp.trigger$depend()
    annot.trigger$depend()
    dimPlotlyOutput(assay.in = input$assay_2, reduc.in = input$reduction_2, group.by = input$grouping_2, annot_panel = input$annot_panel, tmp_annotations = tmp_annotations, low.res = 'yes', data = data)
  })
  
  #-- dimensional reduction plot coloured by feature expression --#
  output$featureplot_2 <- renderPlotly({
    req(data_update())
    req(input$featureplot_2_feature_select)
    tmp.trigger$depend()
    annot.trigger$depend()
    featurePlotlyOutput(assay.in = input$assay_2, reduc.in = input$reduction_2, group.by = input$grouping_2, feature.in = input$featureplot_2_feature_select, low.res = 'yes', data = data)
  })
  
  #-- Find the markers for the 1) selected cells compared to all other cells or
  #-- 2) the markers that define each group. Depending on if cells are selected or not
  observeEvent(input$find.markers, ignoreNULL = FALSE, ignoreInit = FALSE, handlerExpr = {
    output$markers <- DT::renderDataTable({
      req(data_update())
      req(input$grouping_2)
      grouping <- paste0(input$grouping_2,"_meta_data")
      if(!any(names(data[[input$assay_2]]$col.attrs) == grouping)){
        return(NULL)
      }
      if(is.null(assay_matrices) | identical(assay_matrices, list())){
        return(NULL)
      }
      group_assay <- paste(input$grouping_2, input$assay_2, sep = '_')
      if(any(group_assay%in%names(marker_tables)) & is.null(cells)){
        #print('in')
        #t1 <- Sys.time()
        t <- marker_tables[[group_assay]]
        #print(Sys.time()-t1)
      }else{
        #print('not in')
        #t1 <- Sys.time()
        #m <- t(data[[input$assay_2]]$matrix[,])
        #rownames(m) <- data[[input$assay_2]]$row.attrs$features[]
        #colnames(m) <- data[[input$assay_2]]$col.attrs$CellID[]
        y <- data[[input$assay_2]][[paste0('col_attrs/',grouping)]][]
        #print(Sys.time() - t1)
        cell_ids <- colnames(assay_matrices[[input$assay_2]])
        #print(Sys.time() - t1)
        cols <- unlist(lapply(cells[,1], function(x){
          grep(x ,cell_ids, fixed = TRUE)
        }))
        #print(Sys.time() - t1)
        if(!is.null(cells)){
          y[cols] <- 'Selected'
        }
        #print(str(assay_matrices[[input$assay_2]]))
        #print(Sys.time() - t1)
        #print(str(wilcoxauc(X = assay_matrices[[input$assay_2]], y = y)))
        #print(unique(wilcoxauc(X = assay_matrices[[input$assay_2]], y = y)$group))
        #t <- as.data.frame(top_markers(wilcoxauc(X = assay_matrices[[input$assay_2]], y = y), n = 10))[1:10,]
        t <- as.data.frame(wilcoxauc(X = assay_matrices[[input$assay_2]], y = y))
        t <- t[-c(5,6)]
        t$Specificity <- t$pct_in/(t$pct_out+1)
        #print(Sys.time() - t1)
        if(is.null(cells)){
          marker_tables[[group_assay]] <<- t
        }
        #print(Sys.time() - t1)
      }
      
      if(!is.null(cells)){
        #print('sel')
        #t1 <- Sys.time()
        t <- t[which(t$group == "Selected"),]
        #print(str(t))
        #print(Sys.time()-t1)
        #names(t) <- "Markers for Selected Cells"
      }#else{
        #t <- t[mixedsort(t$group),]
        #print(str(t))
        #index <- 1:ncol(t)
        #rank <- grep('rank',colnames(t))
        #index <- index[-rank]
        #t <- t[,c(rank, index)]
      #}
      #print(str(marker_tables))
      t %>% DT::datatable() %>% DT::formatRound(columns=c("avgExpr", "logFC", "pval", "padj", "pct_in", "pct_out", "Specificity"), digits=3)
    })
  })
  
  #-- display selected grouping name --#
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
      for(i in 1:length(data)){
        data[[i]]$add.col.attribute(attributes = new, overwrite = TRUE)
      }
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
        for(i in 1:length(data)){
          data[[i]]$add.col.attribute(attributes = new, overwrite = TRUE)
        } 
        annot.trigger$trigger()
      }
    }
  })
  
  #-- genes to search from for ncbi query --#
  output$gene_query <- renderUI({
    req(data_update())
    req(input$assay_2)
    assay <- input$assay_2
    selectInput(inputId = 'gene_query', label = 'Select a gene of interest to learn more', choices = data[[assay]][['row_attrs/features']][], selected = NA, multiple = FALSE)
  })
  
  #-- get gene summary from ncbi --#
  observeEvent(input$query_ncbi,{
    if(is.null(input$gene_query) | is.na(input$gene_query) | identical(input$gene_query, "") | identical(input$gene_query, " ")){
      showNotification("A valid gene must be entered", type = "error")
      return(NULL)
    }
    id <- fromJSON(file = paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=", input$gene_query, "[sym]+AND+",input$organism,"[orgn]&retmode=json"))$esearchresult$idlist
    #print(id)
    if(identical(id, list())){
      output$gene_summary <- renderText({"Error: No genes matching this query were found on NCBI"})
    }else{
      if(length(id)>1){
        showNotification('Caution: More that one gene matched this query. Only showing first match', type = 'warning')
        id <- id[1]
      }
      output$gene_summary <- renderText({
        results <- fromJSON(file = paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=",id,"&retmode=json"))$result[[2]]
        #print(str(results))
        name <- results$name
        #print(name)
        alias <- results$otheraliases
        #print(alias)
        summary <- results$summary
        #print(summary)
        return(paste0("[id=",id,"][name=",name,"][aliases=", paste(alias, collapse = " "),"][summary=",summary,"]"))
        })
    }
  })
  
  
  
  ###==============// CUSTOM META DATA TAB //==============####
  
  #-- reactive triggers used by this tab --#
  
  # triggered when meta data is selected/deselected from the checkbox group
  meta.trigger <- makeReactiveTrigger()

  #-- global variables used by this tab --#
  
  # orders slected meta data by which order it was selected, rather than alpha-numerically 
  order <- NULL
  
  #-- select input for Assay --#
  output$assay_3 <- renderUI({
    req(data_update())
    selectInput('assay_3', "Select Assay", choices = names(data), selected = ifelse(any(names(data)=="RNA"),yes = "RNA",no = names(data)[1]))
  })
  
  #-- select meta data to combine --#
  output$meta_group_checkbox <- renderUI({
    req(data_update())
    req(input$assay_3)
    choices <- names(data[[input$assay_3]]$col.attrs)[which(grepl("_meta_data$", names(data[[input$assay_3]]$col.attrs)))]
    choices <- sub('_meta_data','',choices)
    checkboxGroupInput(inputId = 'meta_group_checkbox', label = 'Meta Data', choices = choices)
  })
  
  #-- display the order of the selected meta data --#
  output$example_meta_group_name <- renderText({
    req(data_update())
    meta.trigger$depend()
    if(is.null(order) || is.na(order) || length(order) == 0){
      return(NULL)
    }
    paste(sub('_meta_data$','',order),collapse=" + ")
  })
  
  #-- display a sample of what the combined meta data will resemble --#
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
  
  #-- update the 'order' global variable whenever the checkbox group is selected/deselected --#
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
  
  #-- add custom metadata grouping to loom file --#
  observeEvent(input$add_new_meta_group,{
    if(is.null(input$meta_group_checkbox) || length(input$meta_group_checkbox)<2){
      showNotification('At least two meta data groups must be selected!', type = 'error')
    }else if(is.null(input$new_meta_group) || input$new_meta_group == ''){
      showNotification('A name for the meta data must be provided!', type = 'error')
    }else if(paste0(input$new_meta_group,'_meta_data') %in% names(data[[input$assay_3]]$col.attrs)){
      showNotification('The chosen meta data name already exists! Please select a unique name.', type = 'error')
    }else{
      ex <- data[[input$assay_3]]$get.attribute.df(attributes = order)
      ex_list <- list(apply(X = ex, MARGIN = 1, FUN = paste, collapse = '_'))
      if(grepl('_meta_data$', input$new_meta_group)==FALSE){
        names(ex_list) <- paste0(input$new_meta_group,'_meta_data')
      }else{
        names(ex_list) <- input$new_meta_group
      }
      data[[input$assay_3]]$add.col.attribute(attributes = ex_list, overwrite = TRUE)
      update.trigger$trigger()
    }
  })
  
  ###==============// FILE CONVERSION TAB //==============####
  
  shinyFileChoose(input, "object_file", roots = getVolumes(), session = session)
  shinyDirChoose(input, "new_directory", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  output$chosen_object_file <- renderText(parseFilePaths(selection = input$object_file, roots = getVolumes())$datapath)
  
  output$chosen_new_directory <- renderText(parseDirPath(roots = volumes, selection = input$new_directory))
  
  observeEvent(input$object_file, ignoreInit = T, ignoreNULL = T, {
    file_path <- parseFilePaths(selection = input$object_file, roots = getVolumes())$datapath
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
    if(is.null(input$object_file)){
      showNotification('A seurat object must be selected', type = 'error')
      return(NULL)
    }else if(is.null(input$new_directory)){
      showNotification('A directory to save the converted file to must be selected', type = 'error')
      return(NULL)
    }
    dir_path <- parseDirPath(roots = volumes, selection = input$new_directory)

    if(length(list.files(dir_path))>0){
      showNotification('The selected directory already contains files. Please used a new/empty directory.', type = 'error')
      return(NULL)
    }
    
    file_path <- parseFilePaths(selection = input$object_file, roots = getVolumes())$datapath
    
    if(!grepl('\\.rds$', file_path, ignore.case = T)){
      showNotification('The selected file must be of type .rds', type = 'error')
      return(NULL)
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
  
  shinyFileChoose(input, "update_file", roots = getVolumes(), session = session)
  shinyDirChoose(input, "new_directory_2", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  output$chosen_update_file <- renderText(parseFilePaths(selection = input$update_file, roots = getVolumes())$datapath)
  
  output$chosen_new_directory_2 <- renderText(parseDirPath(roots = volumes, selection = input$new_directory_2))
  
  observeEvent(input$ls_convert, ignoreInit = T, {
    if(is.null(input$update_file)){
      showNotification('A seurat object must be selected', type = 'error')
      return(NULL)
    }else if(is.null(input$new_directory_2)){
      showNotification('A directory to save the converted file to must be selected', type = 'error')
      return(NULL)
    }else if(is.null(input$new_file_name) & input$overwrite == "no"){
      showNotification('A file name for the updated objet must be provided', type = 'error')
      return(NULL)
    }
    dir_path <- parseDirPath(roots = volumes, selection = input$new_directory_2)
    path_to_seurat <- parseFilePaths(selection = input$update_file, roots = getVolumes())$datapath
    
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
    error <- loomToSeurat(obj = path_to_seurat, loom = data, dir = dir_path  , file = file)
    removeModal(session = getDefaultReactiveDomain())
    if(error != 1){
      showNotification('Error in file conversion!', type = 'error', closeButton = TRUE, duration = 100)
    }else{
      showNotification('Conversion Complete!', type = 'message', closeButton = F, duration = 15)
    }
  })
  
  
} # server end
