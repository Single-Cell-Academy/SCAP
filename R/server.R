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
library("parallel")

server <- function(input, output, session){
  session$onSessionEnded(stopApp)
  
  options(shiny.maxRequestSize=500*1024^2)
  
  ## Determine folders for ShinyDir button
  #volumes <- c(projects = "/home/joel/SCAP/test_data/projects/) #getVolumes()
  volumes <- c("FTP" = "/ftp",
               Home = fs::path_home())

  ## Source functions
  source("SCAP_functions.R")
  
  # import data
  output$data_used <- renderText({
    if(!is.null(input$assay_1) & !is.na(input$directory[[1]][2])){
      return(paste0("Chosen data: ", input$directory[[1]][2]))
    }else{
      return(paste0("Please select your data..."))
    }
  })
  
  # Old specification of directories
  #shinyDirChoose(input, "directory", roots = volumes, session = session, restrictions = system.file(package = "base"))

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
  shinyDirChoose(input, "directory",
                 roots = volumes,
                 restrictions = system.file(package = "base"))
  
  # connect to .loom files in chosen dir
  loom_data <- reactive({
    req(input$directory)
    path <- parseDirPath(volumes, input$directory)
    if(identical(path, character(0))) return(NULL)
    loom_files <- paste0(path,"/",list.files(path))
    assays <- sub(".loom","",sub(paste0(".*/"),"",loom_files))
    data <- list()
    ## Iterate over all assays and connect to loom objects
    for(i in 1:length(assays)){
      data[[i]] <- tryCatch({
        loomR::connect(loom_files[i], mode = "r+")
      },
        error = function(e){
          showModal(modalDialog(p(paste0("An error occured trying to connect to ", loom_files[i], ". It is possible that someone else is currently analyzing this file. Please wait for them to finish so your annotations do not conflict.")), title = "Error connecting to loom file."), session = getDefaultReactiveDomain())
          return(NULL)
      })
    }
    # if(is.null(data)) return(NULL)
    # if(length(data) != length(assays)) return(NULL)
    # if(length(unlist(lapply(data, function(x){x}))) != length(assays)) return(NULL)
    names(data) <- assays
    return(data)
  })
  
  rvalues <- reactiveValues(tmp_annotations = NULL, cells = NULL, order = NULL)
  
  #-- Make the list of available metadata slots for plotting and comparing a reactive value--#
  loom.trigger <- makeReactiveTrigger() # triggered when loom files are modified, e.g. metadata added to them
  metadata_options <- reactive({
    req(loom_data())
    loom.trigger$depend()
    cols <- names(loom_data()[[1]]$col.attrs)
    grouping_options <- cols[grep('_meta_data$',cols)]
    df <- loom_data()[[1]]$get.attribute.df(attributes = grouping_options)
    pass <- unlist(lapply(grouping_options, function(x){length(unique(df[,x]))<100}))
    grouping_options <- grouping_options[which(pass)]
    return(grouping_options)
    })

  
  marker_tables <- list()
  
  #-- and store normalized count matrix in memory for speed --#
  #=========== SHOULD DECIDE IF DOING THIS IS A GOOD IDEA OR NOT ==========#
  # observeEvent(loom_data(), ignoreInit = TRUE, {
  #   showModal(modalDialog(p("Initializing data for better performance...This may take a little while."), title = paste0("Getting things set up for ", input$directory[[1]][2])), session = getDefaultReactiveDomain())
  #   marker_tables <<- list()
  #   for(x in names(data)){
  #     assay_matrices[[x]] <<- t(Matrix(data = data[[x]]$matrix[,], sparse = TRUE))
  #     rownames(assay_matrices[[x]]) <<- data[[x]]$row.attrs$features[]
  #     colnames(assay_matrices[[x]]) <<- data[[x]]$col.attrs$CellID[]
  #   }
  #   cells <<- NULL
  #   removeModal(session = getDefaultReactiveDomain())
  # })

  ###==============// MAIN TAB //==============####
  
  #-- select input for Assay --#
  output$assay_1 <- renderUI({
    req(loom_data())
    selectInput('assay_1', "Select Assay", choices = names(loom_data()), selected = ifelse(any(names(loom_data())=="RNA"),yes = "RNA",no = names(loom_data())[1]))
  })
  
  #-- select how to group cells --#
  output$grouping_1 <- renderUI({
    req(input$assay_1, metadata_options())
    grouping_options <- metadata_options()
    assay <- input$assay_1
    if(any(grepl(paste0(tolower(assay),"_clusters"),sub("_meta_data","",grouping_options), fixed = TRUE))){
      sel <- paste0(tolower(assay),"_clusters")
    }else{
      sel <- sub("_meta_data","",grouping_options[1])
    }
    selectInput(inputId = 'grouping_1', label = 'Group By', choices = sub("_meta_data$","",metadata_options()), selected = sel, multiple = FALSE)
  })
  
  #-- select the clustering method --#
  output$reduction_1 <- renderUI({
    req(input$assay_1)
    assay <- input$assay_1
    options <- names(loom_data()[[assay]]$col.attrs)
    l.assay <- tolower(assay)
    reductions <- unique(sub("_._reduction$","",options[which(grepl('_reduction$',options))]))
    selectInput(inputId = 'reduction_1', 'Choose Clustering Method', choices = as.list(reductions), selected = ifelse(test = any(reductions==paste0('umap_',l.assay,'_2d')), yes = paste0('umap_',l.assay,'_2d'), no = reductions[1]))
  })
  
  #-- select the featutres for the feature plot --#
  output$featureplot_1_feature.select<- renderUI({
    req(input$assay_1)
    assay <- input$assay_1
    selectInput(inputId = 'featureplot_1_feature.select', label = 'Select a Feature to Visualize on the Feature Plot', choices = loom_data()[[assay]][['row_attrs/features']][], selected = loom_data()[[assay]][['row_attrs/features']][1], multiple = FALSE)
  })
  
  #-- select the featutres for the dot plot --#
  output$dotplot_1_feature_select<- renderUI({
    req(input$assay_1)
    assay <- input$assay_1
    selectInput(inputId = 'dotplot_1_feature_select', label = 'Choose Features to Visualize on Dot Plot', choices = loom_data()[[assay]][['row_attrs/features']][], selected = loom_data()[[assay]][['row_attrs/features']][3], multiple = TRUE)
  })
  
  #-- visualize dotplot as a split dotplot --#
  output$do_split <- renderUI({
    req(input$assay_1)
    radioButtons(inputId = 'do_split', label = 'Split dot plot', choices = c('yes', 'no'), selected = 'no', inline = TRUE)
  })
  
  #-- select how to split the dot plot --#
  output$split_by <- renderUI({
    req(input$assay_1, metadata_options())
    choices <- unlist(lapply(metadata_options(), function(x){
      if(length(unique(loom_data()[[input$assay_1]]$col.attrs[[x]][]))>1){
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
    req(input$grouping_1, input$reduction_1)
    dimPlotlyOutput(assay.in = input$assay_1, reduc.in = input$reduction_1, group.by = input$grouping_1, annot_panel = "", low.res = 'yes', data = loom_data())
  })
  
  #-- dimensional reduction plot coloured by feature expression --#
  output$featureplot_1 <- renderPlotly({
    req(input$featureplot_1_feature.select)
    featurePlotlyOutput(assay.in = input$assay_1, reduc.in = input$reduction_1, group.by = input$grouping_1, feature.in = input$featureplot_1_feature.select, low.res = 'yes', data = loom_data())
  })
  
  #-- dot plot for selected feature expression --#
  dot_plot <- reactive({
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
      p <- split_dot_plot(data = loom_data()[[assay]], features = input$dotplot_1_feature_select, group.by = group.by, split.by = input$split_by)
      if(is.null(p)) return(NULL)
      if(identical(class(p),"shiny.tag")) return(NULL)
      return(p)
    }else if(input$do_split == 'no'){
      p <- dotPlot(data = loom_data()[[assay]], assay = assay, features = input$dotplot_1_feature_select, group.by = group.by)
      if(is.null(p)) return(NULL)
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
    req(loom_data())
    selectInput('assay_2', "Select Assay", choices = names(loom_data()), selected = ifelse(any(names(loom_data())=="RNA"),yes = "RNA",no = names(loom_data())[1]))
  })
  
  #-- initialize/reset tmp_annotations --#
  observeEvent(input$assay_2,{
    req(input$assay_2)
    rvalues$tmp_annotations <- rep("unlabled", loom_data()[[input$assay_2]]$shape[2])
    names(rvalues$tmp_annotations) <- loom_data()[[input$assay_2]]$col.attrs$CellID[]
  })
  
  #-- select how to group cells --#
  output$grouping_2 <- renderUI({
    req(input$assay_2, metadata_options())
    grouping_options <- metadata_options()
    assay <- input$assay_2
    if(any(grepl(paste0(tolower(assay),"_clusters"),sub("_meta_data","",grouping_options), fixed = TRUE))){
      sel <- paste0(tolower(assay),"_clusters")
    }else{
      sel <- sub("_meta_data","",grouping_options[1])
    }
    selectInput(inputId = 'grouping_2', label = 'Group By', choices = sub("_meta_data$","",metadata_options()), selected = sel, multiple = FALSE)
  })
  
  #-- select the clustering method --#
  output$reduction_2 <- renderUI({
    req(input$assay_2)
    assay <- input$assay_2
    options <- names(loom_data()[[assay]]$col.attrs)[grep("_2d",names(loom_data()[[assay]]$col.attrs))]
    l.assay <- tolower(assay)
    reductions <- unique(sub("_._reduction$","",options[which(grepl('_reduction$',options))]))
    selectInput(inputId = 'reduction_2', 'Choose Clustering Method', choices = as.list(reductions), selected = ifelse(test = any(reductions==paste0('umap_',l.assay,'_2d')), yes = paste0('umap_',l.assay,'_2d'), no = reductions[1]))
  })
  
  #-- select the featutres for the feature plot --#
  output$featureplot_2_feature_select<- renderUI({
    req(input$assay_2)
    assay <- input$assay_2
    choices <- loom_data()[[assay]][['row_attrs/features']][]
    selectInput(inputId = 'featureplot_2_feature_select', label = 'Select a Feature to Visualize on the Feature Plot', choices = choices, selected = choices[1], multiple = FALSE)
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
    req(input$assay_2, metadata_options())
    dimPlotlyOutput(assay.in = input$assay_2, reduc.in = input$reduction_2, group.by = input$grouping_2, annot_panel = input$annot_panel, tmp_annotations = rvalues$tmp_annotations, low.res = 'yes', data = loom_data())
  })
  
  #-- dimensional reduction plot coloured by feature expression --#
  output$featureplot_2 <- renderPlotly({
    req(input$featureplot_2_feature_select)
    featurePlotlyOutput(assay.in = input$assay_2, reduc.in = input$reduction_2, group.by = input$grouping_2, feature.in = input$featureplot_2_feature_select, low.res = 'yes', data = loom_data())
  })
  
  #-- Find the markers for the 1) selected cells compared to all other cells or
  #-- 2) the markers that define each group. Depending on if cells are selected or not
  observeEvent(input$find.markers, ignoreNULL = FALSE, ignoreInit = FALSE, handlerExpr = {
    output$markers <- DT::renderDataTable({
      req(input$assay_2, metadata_options())
      cells <- isolate(rvalues$cells)
      grouping <- paste0(input$grouping_2,"_meta_data")
      if(!any(names(loom_data()[[input$assay_2]]$col.attrs) == grouping)){
        return(NULL)
      }
      group_assay <- paste(input$grouping_2, input$assay_2, sep = '_')
      if(any(group_assay%in%names(marker_tables)) & is.null(cells)){
        t <- marker_tables[[group_assay]]
      }else{
        y <- loom_data()[[input$assay_2]][[paste0('col_attrs/',grouping)]][]
        cell_ids <- loom_data()[[input$assay_2]]$col.attrs$CellID[]
        cols <- unlist(lapply(cells[,1], function(x){
          grep(x ,cell_ids, fixed = TRUE)
        }))
        if(!is.null(cells)){
          y[cols] <- 'Selected'
        }
        #t <- as.data.frame(top_markers(wilcoxauc(X = assay_matrices[[input$assay_2]], y = y), n = 10))[1:10,]
        exp_mat <- Matrix::t(Matrix::Matrix(loom_data()[[input$assay_2]]$matrix[,], sparse = TRUE))
        rownames(exp_mat) <- loom_data()[[input$assay_2]]$row.attrs$features[]
        colnames(exp_mat) <- loom_data()[[input$assay_2]]$col.attrs$CellID[]
        t <- as.data.frame(wilcoxauc(X = exp_mat, y = y))
        t <- t[-c(5,6)]
        t$Specificity <- t$pct_in/(t$pct_out+1)
        if(is.null(cells)){
          marker_tables[[group_assay]] <<- t
        }
      }
      if(!is.null(cells)){
        t <- t[which(t$group == "Selected"),]
      }
      return(t %>% arrange(desc(Specificity)) %>% DT::datatable(filter = 'top') %>% DT::formatRound(columns=c("avgExpr", "logFC", "pval", "padj", "pct_in", "pct_out", "Specificity"), digits=3))
    })
  })
  
  #-- display selected grouping name --#
  output$annotation_title <- renderText({
    req(input$grouping_2, metadata_options())
    return(paste0(input$grouping_2))
  })
  
  #-- display text boxes for user to enter new IDs for the meta data --#
  output$annotations <- renderUI({
    req(input$assay_2, metadata_options())
    group.by <- paste0(input$grouping_2,'_meta_data')
    if(!group.by%in%names(loom_data()[[input$assay_2]]$col.attrs)){
      return(NULL)
    }
    names <- unique(loom_data()[[input$assay_2]][[paste0('col_attrs/',group.by)]][])
    names <- mixedsort(names)
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

  #-- save user defined annotations --#
  observeEvent(input$set_cell_types, ignoreNULL = FALSE, ignoreInit = FALSE, handlerExpr = {
    req(metadata_options())
    group.by <- paste0(input$grouping_2,'_meta_data')
    id <- "my_annot"
    if(is.null(input$new_scheme_name) || input$new_scheme_name == "" || input$new_scheme_name == " "){
      showNotification("Please name the annotation scheme", type = 'warning')
    }else if(input$annot_panel == 'cell_annotation_cluster'){
      names <- loom_data()[[input$assay_2]][[paste0('col_attrs/',group.by)]][]
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
      for(i in 1:length(loom_data())){
        loom_data()[[i]]$add.col.attribute(attributes = new, overwrite = TRUE)
      }
    }else{
      if(all(rvalues$tmp_annotations=='unlabled')){
        showNotification("No annotations have been added", type = 'warning')
      }else if(is.null(input$new_scheme_name) | input$new_scheme_name==""){
        showNotification("Please enter an name for the annotation scheme", type = 'warning')
      }else{
        new <- list(rvalues$tmp_annotations)
        names(new) <- input$new_scheme_name
        if(!grepl("_meta_data$",names(new))){
          names(new) <- paste0(names(new),"_meta_data")
        }
        for(i in 1:length(loom_data())){
          loom_data()[[i]]$add.col.attribute(attributes = new, overwrite = TRUE)
        } 
      }
    }
    loom.trigger$trigger()
  })
  
  #-- genes to search from for ncbi query --#
  output$gene_query <- renderUI({
    req(input$assay_2)
    assay <- input$assay_2
    selectInput(inputId = 'gene_query', label = 'Select a gene of interest to learn more', choices = loom_data()[[assay]][['row_attrs/features']][], selected = NA, multiple = FALSE)
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
    
  #-- select input for Assay --#
  output$assay_3 <- renderUI({
    req(loom_data())
    selectInput('assay_3', "Select Assay", choices = names(loom_data()), selected = ifelse(any(names(loom_data())=="RNA"),yes = "RNA",no = names(loom_data())[1]))
  })
  
  #-- select meta data to combine --#
  output$meta_group_checkbox <- renderUI({
    req(input$assay_3, loom_data())
    choices <- metadata_options()
    choices <- sub('_meta_data','',choices)
    checkboxGroupInput(inputId = 'meta_group_checkbox', label = 'Meta Data', choices = choices)
  })
  
  #-- display the order of the selected meta data --#
  output$example_meta_group_name <- renderText({
    req(loom_data())
    if(is.null(rvalues$order) || is.na(rvalues$order) || length(rvalues$order) == 0){
      return(NULL)
    }
    paste(sub('_meta_data$','',rvalues$order),collapse=" + ")
  })
  
  #-- display a sample of what the combined meta data will resemble --#
  output$example_meta_group <- renderTable({
    req(input$assay_3)
    if(is.null(rvalues$order) || is.na(rvalues$order) || length(rvalues$order) == 0){
      return(NULL)
    }
    ex <- loom_data()[[input$assay_3]]$get.attribute.df(attributes = rvalues$order)
    ex_df <- as.data.frame(apply(X = ex, MARGIN = 1, FUN = paste, collapse = '_'), stringsAsFactors = FALSE)[sample(x = 1:nrow(ex),size = 20,replace = FALSE),,drop=F]
    colnames(ex_df) <- 'Examples of chosen meta data combination'
    return(ex_df)
  }, hover = TRUE, width = '100%',align = 'c')
  
  #-- update the 'order' global variable whenever the checkbox group is selected/deselected --#
  observeEvent(input$meta_group_checkbox, ignoreNULL = F,{
    req(loom_data())
    if(!is.null(input$meta_group_checkbox)){
      new <- paste0(input$meta_group_checkbox,'_meta_data')
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
    }else if(paste0(input$new_meta_group,'_meta_data') %in% names(loom_data()[[input$assay_3]]$col.attrs)){
      showNotification('The chosen meta data name already exists! Please select a unique name.', type = 'error')
    }else{
      ex <- loom_data()[[input$assay_3]]$get.attribute.df(attributes = rvalues$order)
      ex_list <- list(apply(X = ex, MARGIN = 1, FUN = paste, collapse = '_'))
      if(grepl('_meta_data$', input$new_meta_group)==FALSE){
        names(ex_list) <- paste0(input$new_meta_group,'_meta_data')
      }else{
        names(ex_list) <- input$new_meta_group
      }
      showModal(modalDialog(p(paste0("Adding ", input$new_meta_group, " to meta data...")), title = "This window will close after once complete"), session = getDefaultReactiveDomain())
      Sys.sleep(2)
      for(i in 1:length(loom_data())){
        loom_data()[[i]]$add.col.attribute(attributes = ex_list, overwrite = TRUE)
      }
      loom.trigger$trigger()
      Sys.sleep(2)
      removeModal(session = getDefaultReactiveDomain())
    }
  })
  
  ###==============// FILE CONVERSION TAB //==============####
  
  shinyFileChoose(input, "object_file", roots = volumes, session = session)
  shinyDirChoose(input, "new_directory", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  output$chosen_object_file <- renderText(parseFilePaths(selection = input$object_file, roots = volumes)$datapath)
  
  output$chosen_new_directory <- renderText(parseDirPath(roots = volumes, selection = input$new_directory))
  
  observeEvent(input$object_file, ignoreInit = T, ignoreNULL = T, {
    print(input$object_file)
    message(input$object_file)
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
    
    file_path <- parseFilePaths(selection = input$object_file, roots = volumes)$datapath
    
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
  
  shinyFileChoose(input, "update_file", roots = volumes, session = session)
  shinyDirChoose(input, "new_directory_2", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
  output$chosen_update_file <- renderText(parseFilePaths(selection = input$update_file, roots = volumes)$datapath)
  
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
  

  ####
  #-- Functions for annotation comparison--#
  ####

  ## Make meta data annotations a reactive object instead of performing all of the calles 3 times! 

  #-- get the first list of potential annotations from loom--#
  output$comp_anno_list1 <- renderUI({
    req(loom_data(), metadata_options())

    grouping_options <- metadata_options()
    
    selectInput(
      inputId = 'comp_anno_1', 
      label = 'Select the annotation you want to compare!', 
      choices = sub("_meta_data$","",metadata_options()), 
      multiple = FALSE)
  })
  
  #-- second list of annotations to compare against--#
  output$comp_anno_list2 <- renderUI({
    req(loom_data(), input$comp_anno_1)

    grouping_options <- metadata_options()

    grouping_options <- grouping_options[!grepl(paste(input$comp_anno_1,"_meta_data",sep=""),grouping_options)]
    
    selectInput(
      inputId = 'comp_anno_2', 
      label = 'Select the annotation to compare against!', 
      choices = sub("_meta_data$","",grouping_options), 
      multiple = FALSE)  })
  


  sankey_comp <- eventReactive(input$compare_annos, {
    req(input$comp_anno_1, input$comp_anno_2)

    annos_to_compare <- loom_data()[[input$assay_1]]$get.attribute.df(
      attributes=c(
        "CellID",
        paste(input$comp_anno_1,"_meta_data",sep=""),
        paste(input$comp_anno_2,"_meta_data",sep=""))
    )

    colnames(annos_to_compare) <- gsub("_meta_data","",colnames(annos_to_compare))
    
    annos_to_compare_stats <- annos_to_compare %>% 
      group_by(get(input$comp_anno_1),get(input$comp_anno_2)) %>%
      tally() %>%
      ungroup()
    
    colnames(annos_to_compare_stats) <- c("anno1","anno2","n")
    
    annos_to_compare_stats <- annos_to_compare_stats %>%
      mutate("anno1" = paste(anno1,input$comp_anno_1,sep="_")) %>%
      mutate("anno2" = paste(anno2,input$comp_anno_2,sep="_"))
    
    joined_annos <- c(annos_to_compare_stats$anno1,annos_to_compare_stats$anno2)
    joined_annos <- unique(joined_annos)
    
    annos_to_compare_stats$IDsource=match(annos_to_compare_stats$anno1, joined_annos)-1 
    annos_to_compare_stats$IDtarget=match(annos_to_compare_stats$anno2, joined_annos)-1

    return(annos_to_compare_stats)
    })
  


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
        thickness = 20,
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
    )



    }
  )
  
} # server end
