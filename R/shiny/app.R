#################### load required packages ###############
pkgs <- c('shiny','Seurat','ggplot2','ggthemes','BiocManager','plotly','cowplot', 'dplyr','devtools','shinyjqui','shinythemes')
lapply(X = pkgs, FUN = function(x){
  if(!require(x, character.only = TRUE))
    install.packages(x, repos = 'https://cloud.r-project.org/')
})

lapply(X = pkgs, FUN = function(x){
  if(!require(x, character.only = TRUE))
    library(x)
})
if(!require('shinycssloaders'))
  devtools::install_github('https://github.com/andrewsali/shinycssloaders')
library(shinycssloaders)
if(!require('presto'))
  devtools::install_github('https://github.com/immunogenomics/presto.git')
library(presto)
library(shinycssloaders)
if(!require('loomR'))
  BiocManager::install('loomR')
library(loomR)



#---------------------------------------------------------#
#                   User Interface                        #
#_________________________________________________________#

ui <- navbarPage(fluid = TRUE,
  theme = shinytheme('cosmo'),
  title = "Single Cell Analysis Portal",
  tabPanel(
    "Main",
    sidebarPanel(
      "",
      uiOutput('assay.1') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '100px'),
      uiOutput('reduction.1') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '100px'),
      uiOutput('grouping.1') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '100px'),
      uiOutput('featureplot.1.feature.select') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '100px'),
      uiOutput('dotplot.1.feature.select') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '100px')
    )#,
    # mainPanel(
    #   "",
    #   fluidRow(style='padding:100px;',
    #     column(
    #       width = 6,
    #       plotlyOutput('dimplot.1', height = '100%', width = '100%')%>% withSpinner(color="black")
    #     ),
    #     column(
    #       width = 6,
    #       plotlyOutput('featureplot.1', height = '100%', width = '100%')%>% withSpinner(color="black")
    #     ),
    #     div(style = "height:10px;")
    #   ),
    #   fluidRow(
    #     column(
    #       width = 12, 
    #       offset = 0,
    #       plotOutput('dotplot.1', height = '500px', width = '100%')%>% withSpinner(color="black")
    #     )
    #   )
    # )
  )#,
  # tabPanel(
  #   "Cell Annotation",
  #   sidebarPanel(
  #     uiOutput('assay.2'),
  #     uiOutput('reduc.2'),
  #     uiOutput('feature.select3'),
  #     tabsetPanel(
  #       id = 'annot_panel',
  #       tabPanel(
  #         value = 'cluster',
  #         "Annotate Clusters",
  #         uiOutput('group.by.2'),
  #         uiOutput('annots'),
  #         textInput('new_name', 'Name your annotation scheme', placeholder = "Enter scheme name here..."),
  #         actionButton("setCellsTypes", "Rename Clusters/Cells")
  #       ),
  #       tabPanel(
  #         value = 'custom',
  #         "Custom Annotations",
  #         textInput('custom_name','Name Selected Cluster', placeholder = "Enter cluster name here..."),
  #         actionButton("add_to_tmp", "Add Annotation"),
  #         textInput('custom_scheme','Name Custom Scheme', placeholder = "Enter scheme name here..."),
  #         actionButton("save_custom_scheme", "Save Custom Annotations")
  #       )
  #     ),
  #     h4("Selected Cells"),
  #     tableOutput('cells'),
  #     tags$head(tags$style("#cells{overflow-y:scroll; max-height: 200px; background: ghostwhite;}"))
  #   ),
  #   mainPanel(
  #     fluidRow(style='padding:100px;',
  #       column(width = 6,
  #              plotlyOutput('annot.select.cluster', height = '100%', width = '100%')%>% withSpinner(color="black")
  #              ),
  #       column(width = 6,
  #              plotlyOutput('annot.select.feature', height = '100%', width = '100%')%>% withSpinner(color="black")
  #              )
  #     ),
  #     fluidRow(style='padding:50px;',
  #       actionButton('find.markers', "Find Markers for Selected Cells"),
  #       h4('To clear selection... Double click on either plot and click button')
  #       ),
  #     fluidRow(
  #       tableOutput('markers') %>% withSpinner(color="black")
  #     )
  #   )
  # ),
  # tabPanel(
  #   "Create Custom Meta Data Groupings",
  #   sidebarPanel(
  #     uiOutput(outputId = 'assay.3'),
  #     uiOutput(outputId = 'meta_group_checkbox'),
  #     textInput(inputId = 'new_meta_group', placeholder = 'Enter meta data grouping name here...',label = NULL),
  #     actionButton(inputId = 'add_new_meta_group', label = "Add")
  #   ),
  #   mainPanel(
  #     h1(textOutput(outputId = 'example_meta_group_name')),
  #     textOutput(outputId = 'example_meta_group')
  #   )
  # )
)
#jqui_resizable(jqui_draggable(
#%>% withSpinner(color="#0dc5c1")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


#---------------------------------------------------------#
#                       SERVER                            #
#_________________________________________________________#

server <- function(input, output, session){

  ###### set working dir and source required functions #####
  options(shiny.maxRequestSize=500*1024^2)
  setwd("/Users/jsjoyal/Desktop/SCAP/")
  source("./R/SCAP_functions.R")
  
  ###### import data ######
  files <- list.files("./test_data/test_project_4/")
  assays <- assays <- sub(".loom","",files)
  data <- list()
  for(i in 1:length(assays)){
    data[[i]] <- connect(paste0("./test_data/test_project_4/",files[i]), mode = "r+")
  }
  names(data) <- assays
  
  print(data)
  
  ## Golbal Variables ##
  grouping_options_old <- NULL
  
  ####################### MAIN TAB ########################
  #                                                       #
  ####################### MAIN TAB ########################
  

  ########## Select Input for Assay ##############
  output$assay.1 <- renderUI({selectInput('assay.1', "Select Assay", choices = names(data), selected = ifelse(any(names(data)=="RNA"),yes = "RNA",no = names(data)[1]))})
  
  output$reduction.1 <- renderUI({
    req(input$assay.1)
    assay <- input$assay.1
    #print(assay)
    #reduction
    options <- names(data[[assay]]$col.attrs)
    l.assay <- tolower(assay)
    reductions <- unique(sub("_._reduction$","",options[which(grepl('_reduction$',options)&grepl(l.assay,options))]))
    selectInput(inputId = 'reduction.1', 'Select Reduction', choices = as.list(reductions), selected = ifelse(test = any(reductions==paste0('tsne_',l.assay,'_2d')), yes = paste0('tsne_',l.assay,'_2d'), no = reductions[1]))
  })

  output$grouping.1 <- renderUI({
    req(input$assay.1)
    assay <- input$assay.1
    #grouping
    cols <- names(data[[assay]]$col.attrs)
    grouping_options <- cols[grep('_meta_data$',cols)]
    if(!is.null(grouping_options_old) & (length(grouping_options)!=length(grouping_options_old))){
      sel <- sub("_meta_data$", "", grouping_options[!grouping_options%in%grouping_options_old])
    }else{
      sel <- sub("_meta_data$", "", grouping_options)
    }
    grouping_options_old <<- grouping_options

    selectInput(inputId = 'grouping.1', label = 'Choose Annotations', choices = sub("_meta_data$","",grouping_options), selected = sel, multiple = FALSE)
  })
  # 
  # output$featureplot.1.feature.select<- renderUI({
  #   req(input$assay.1)
  #   assay <- input$assay.1
  #   selectInput(inputId = 'featureplot.1.feature.select', label = 'Choose Features for Feature Plot', choices = data[[assay]][['row_attrs/features']][], selected = data[[assay]][['row_attrs/features']][1], multiple = FALSE)
  # })
  # output$dotplot.1.feature.select<- renderUI({
  #   req(input$assay.1)
  #   assay <- input$assay.1
  #   selectInput(inputId = 'dotplot.1.feature.select', label = 'Choose Features for Dot Plot', choices = data[[assay]][['row_attrs/features']][], selected = data[[assay]][['row_attrs/features']][3], multiple = TRUE)
  # })
  # 
  ############### Dim Plot ################
  output$dimplot.1 <- renderPlotly({
    req(input$grouping.1)
    dimPlotlyOutput(assay.in = input$assay.1, reduc.in = input$reduction.1, group.by = input$grouping.1, data = data)
  })

  # ########## Feature Plot #####################
  # output$featureplot.1 <- renderPlotly({
  #   req(input$featureplot.1.feature.select)
  #   featurePlotlyOutput(assay.in = input$assay.1, reduc.in = input$reduction.1, group.by = input$grouping.1, feature.in = input$featureplot.1.feature.select, data = data)
  # })
  # 
  # ######### Dot Plot #########
  # output$dotplot.1 <- renderPlot({
  #   req(input$dotplot.1.feature.select)
  # 
  #   group.by <- input$grouping.1
  #   assay <- input$assay.1
  # 
  #   ax.x <- list(
  #     title = "Features",
  #     zeroline = FALSE,
  #     showline = TRUE,
  #     showticklabels = TRUE,
  #     showgrid = FALSE,
  #     mirror=TRUE,
  #     ticks='outside',
  #     tickangle = -45
  #   )
  #   ax.y <- list(
  #     title = "Identity",
  #     zeroline = FALSE,
  #     showline = TRUE,
  #     showticklabels = TRUE,
  #     showgrid = FALSE,
  #     mirror=TRUE,
  #     ticks='outside'
  #   )
  #   
  #   p <- dotPlot(data = data[[assay]], assay = assay, features = input$dotplot.1.feature.select, group.by = group.by)
  #   
  #   if(is.null(p)) return(NULL)
  #   if(identical(class(p),"shiny.tag")) return(NULL)
  #   
  #   p <- p + guides(color = guide_colourbar(order = 1, title = "Average Expression"), size = guide_legend(order = 2, title = "Percent Expressed")) + theme_few()
  #   
  #   return(p)
  # })
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  ######### Select input for reduction ##########

  # output$reduction.1 <- renderUI({
  #   req(input$assay.1)
  #   assay <- input$assay.1
  #   id <- 'reduction.1'
  #   options <- c("pca", "tsne", "umap", "diff")
  #   l.assay <- tolower(assay)
  #   options <- paste(c("pca", "tsne", "umap", "diff"), l.assay,sep = "_")
  #   reductions <- unique(sub("_.$","",names(data[[assay]]$col.attrs)[grepl(paste(options,collapse="|"),names(data[[assay]]$col.attrs))]))
  #   print(reductions)
  #   selectInput(
  #     inputId = id,
  #     "Select Reduction",
  #     choices = as.list(reductions),
  #     selected = reductions[1]
  #   )
  # })
  # 
  # ########## Select group.by ###########
  # output$grouping.1 <- renderUI({
  #   req(input$assay.1)
  #   assay <- input$assay.1
  #   id <- 'grouping.1'
  #   cell_type_options <- names(data[[]]$col.attrs)[grep('cell_type|cell types|cell type|cell_types',names(data[[assay]]$col.attrs))]
  #   if(!is.null(cell_type_options_old) & (length(cell_type_options)!=length(cell_type_options_old))){
  #     sel <- cell_type_options[!cell_type_options%in%cell_type_options_old]
  #   }else{
  #     sel <- 'cluster'
  #   }
  #   cell_type_options_old <<- cell_type_options
  #   selectInput(inputId = id, label = 'Choose Annotations', choices = c('cluster', cell_type_options), selected = sel, multiple = FALSE)
  # })


  ########## Select features for dotplot ###########
  # observe({
  #   req(input$assay.1)
  #   assay <- input$assay.1
  #   id <- 'dotplot.1.feature.select'
  #   output[[id]] <- renderUI({
  #     selectInput(inputId = id, label = 'Choose Features for Dot Plot', choices = data[[assay]][['row_attrs/features']][], selected = data[[assay]][['row_attrs/features']][1:3], multiple = TRUE)
  #   })
  # })
  
  ########## Select features for feature plot ###########
  # observe({
  #   req(input$grouping.1)
  #   assay <- input$assay.1
  #   id <- 'featureplot.1.feature.select'
  #   output[[id]] <- renderUI({
  #     selectInput(inputId = id, label = 'Choose a Feature for Feature Plot', choices = data[[assay]][['row_attrs/features']][], selected = data[[assay]][['row_attrs/features']][1], multiple = FALSE)
  #   })
  # })
  

  ####################### ANNOTATION TAB ########################
  #                                                             #
  ####################### ANNOTATION TAB ########################
  
  # tmp_annotations <- NULL
  # 
  # observeEvent(input$assay.2,{
  #   tmp.trigger$depend()
  #   tmp_annotations <<- rep("unlabled",data[[input$assay.2]]$shape[2])
  #   names(tmp_annotations) <<- data[[input$assay.2]]$col.attrs$CellID[]
  # })
  # 
  # tmp.trigger <- makeReactiveTrigger()
  # annot.trigger <- makeReactiveTrigger()
  # 
  # output$annot.select.cluster <- renderPlotly({
  #   tmp.trigger$depend()
  #   annot.trigger$depend()
  #   req(input$reduc.2)
  #   req(input$group.by.2)
  #   req(input$assay.2)
  #   
  #   reduc <- tolower(paste0(input$reduc.2,"_rna_2d_",1:2))
  #   
  #   if(input$group.by.2 == "cluster"){
  #     plot.data <- data[[input$assay.2]]$get.attribute.df(attributes=c(reduc,tolower(paste0(input$assay.2,"_","clusters")),'percent.mt'))
  #   }
  #   else{
  #     plot.data <- data[[input$assay.2]]$get.attribute.df(attributes=c(reduc,input$group.by.2,'percent.mt'))
  #   }
  #   key <- reduc_key(key = toupper(input$reduc.2))
  #   
  #   ax.x <- list(
  #     title = paste0(key,"_1"),
  #     zeroline = FALSE,
  #     showline = TRUE,
  #     showticklabels = FALSE,
  #     showgrid = FALSE,
  #     mirror=TRUE,
  #     ticks='none'
  #   )
  #   ax.y <- list(
  #     title = paste0(key,"_2"),
  #     zeroline = FALSE,
  #     showline = TRUE,
  #     showticklabels = FALSE,
  #     showgrid = FALSE,
  #     mirror=TRUE,
  #     ticks='none'
  #   )
  #   
  #   ax.z <- list(
  #     title = paste0(key,"_3"),
  #     zeroline = FALSE,
  #     showline = TRUE,
  #     showticklabels = FALSE,
  #     showgrid = FALSE,
  #     mirror=TRUE,
  #     ticks='none'
  #   )
  #   
  #   #print(head(plot.data))
  #   
  #   if(input$annot_panel == 'cluster'){
  #     col <- plot.data[,3]
  #   }else{
  #     col <- tmp_annotations
  #   }
  #   p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2], key = ~rownames(plot.data),
  #                type = 'scatter', mode = 'markers', 
  #                color = col, text =  ~paste0(
  #                  key,"_1: ", format(plot.data[,1],digits=3),"\n",
  #                  "</br>",key,"_2: ", format(plot.data[,2],digits=3), "\n",
  #                  "</br>percent.mt: ", format(plot.data[,4],digits=3), "%"), 
  #                hovertemplate = paste0('<b>%{text}</b>')
  #   ) %>% layout(title = ifelse(test = input$annot_panel == 'cluster', yes = paste0(input$assay.2, " data coloured by ", input$group.by.2), no = paste0(input$assay.2, " data")),xaxis = ax.x, yaxis = ax.y, dragmode = "select")
  # 
  #   p
  # })
  # output$annot.select.feature <- renderPlotly({
  #   tmp.trigger$depend()
  #   annot.trigger$depend()
  #   req(input$feature.select3)
  #   req(input$reduc.2)
  #   req(input$assay.2)
  #   data.features <- as.data.frame(data[[input$assay.2]][['matrix']][,which(data[[input$assay.2]]$row.attrs$features[]%in%input$feature.select3)])
  #   if(nrow(data.features)==0|ncol(data.features)==0) return(NULL)
  #   if(input$group.by.2 == 'cluster'){
  #     data.annot <- data[[input$assay.2]]$get.attribute.df(attributes=tolower(paste0(input$assay.2,"_clusters")))
  #   }else{
  #     data.annot <- data[[input$assay.2]]$get.attribute.df(attributes=input$group.by.2)
  #   }
  #   rownames(data.features) <- data[[input$assay.2]]$col.attrs$CellID[]
  #   colnames(data.features) <- input$feature.select3
  # 
  #   
  #   dims <- tolower(paste0(input$reduc.2,"_",input$assay.2,"_2d_",1:2))
  # 
  #   plot.data <- data[[input$assay.2]]$get.attribute.df(attributes=dims)
  #   plot.data <- cbind(plot.data,data.annot)
  #   plot.data <- cbind(plot.data,data.features)
  # 
  #   #plot.data[,4] <- scale(plot.data[,4])
  # 
  #   key <- reduc_key(key = toupper(input$reduc.2))
  # 
  #   ax.x <- list(
  #     title = paste0(key,"_1"),
  #     zeroline = FALSE,
  #     showline = TRUE,
  #     showticklabels = FALSE,
  #     showgrid = FALSE,
  #     mirror=TRUE,
  #     ticks='none'
  #   )
  #   ax.y <- list(
  #     title = paste0(key,"_2"),
  #     zeroline = FALSE,
  #     showline = TRUE,
  #     showticklabels = FALSE,
  #     showgrid = FALSE,
  #     mirror=TRUE,
  #     ticks='none'
  #   )
  #   ax.z <- list(
  #     title = paste0(key,"_3"),
  #     zeroline = FALSE,
  #     showline = TRUE,
  #     showticklabels = FALSE,
  #     showgrid = FALSE,
  #     mirror=TRUE,
  #     ticks='none'
  #   )
  #   plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2],
  #           type = 'scatter', mode = 'markers', key = ~rownames(plot.data),
  #           color = plot.data[,4],text =  ~paste0(
  #             key,"_1: ", format(plot.data[,1],digits=3),"\n",
  #             "</br>",key,"_2: ", format(plot.data[,2],digits=3), "\n",
  #             ifelse(test = input$annot_panel == 'cluster', yes = paste0("</br>",input$group.by.2,": ", plot.data[,3]), no = "")),
  #           hovertemplate = paste0('<b>%{text}</b>',
  #                                  '<extra></extra>')
  #   ) %>% layout(title = input$feature.select3 ,xaxis = ax.x, yaxis = ax.y, dragmode = "select")
  # })
  # 
  # cells <- NULL
  # output$cells <- renderTable({
  #   d <- as.data.frame(event_data("plotly_selected")$key, stringsAsFactors = FALSE)
  #   #print(head(d))
  #   #print(str(d))
  #   if (is.null(d) | ncol(d) == 0 | nrow(d) == 0){
  #     cells <<- NULL
  #     "Click and drag on either plot to select cells for marker identification (i.e., select/lasso) (double-click to clear selection)"
  #   }else{
  #     if(ncol(d)>1){
  #       d <- as.data.frame(as.character(d[1,]),stringsAsFactors = F)
  #     }
  #     colnames(d) <- ""
  #     #print(d)
  #     cells <<- d
  #     d
  #   }
  # })
  # 
  # observeEvent(input$find.markers, ignoreNULL = FALSE, ignoreInit = FALSE, handlerExpr = {
  #   output$markers <- renderTable({
  #     req(input$assay.2)
  #     req(input$group.by.2)
  #     #if(is.null(input$find.markers)) NULL
  #     m <- t(data[[input$assay.2]]$matrix[,])
  #     rownames(m) <- data[[input$assay.2]]$row.attrs$features[]
  #     colnames(m) <- data[[input$assay.2]]$col.attrs$CellID[]
  #     if(input$group.by.2=='cluster'){
  #       y <- data[[input$assay.2]][[paste0("col_attrs/",tolower(input$assay.2),"_clusters")]][]
  #     }else{
  #       y <- data[[input$assay.2]][[paste0('col_attrs/',input$group.by.2)]][]
  #     }
  #     cols <- unlist(lapply(cells[,1], function(x){
  #       grep(x,data[[input$assay.2]]$col.attrs$CellID[])
  #     }))
  #     #print(cells)
  #     #print(cols)
  #     #print(head(cols))
  #     if(!is.null(cells))
  #       y[cols] <- 'Selected'
  #     #print(head(y))
  #     t <- as.data.frame(top_markers(wilcoxauc(X = m, y = y), n = 10))[1:10,]
  #     if(!is.null(cells)){
  #       # t1 <- t[,-which(colnames(t)=="Selected")]
  #       # t <- cbind(t1,t[,which(colnames(t)=="Selected")])
  #       # colnames(t)[ncol(t)] <- "Selected"
  #       t <- as.data.frame(t[,which(colnames(t)=="Selected")],stringsAsFactors=FALSE)
  #       names(t) <- "Markers for Selected Cells"
  #     }
  #     t
  #   }, hover = TRUE, width = '100%',align = 'c')
  # })
  # 
  # observeEvent(input$setCellsTypes, ignoreNULL = FALSE, ignoreInit = FALSE, handlerExpr = {
  #   req(input$assay.2)
  #   req(input$group.by.2)
  #   if(is.null(input$new_name) || input$new_name == "" || input$new_name == " "){
  #     showNotification("Please name the annotation scheme", type = 'warning')
  #   }else{
  #     if(input$group.by.2 == "cluster"){
  #       names <- data[[input$assay.2]][[paste0("col_attrs/",tolower(input$assay.2),"_clusters")]][]
  #       new <- names
  #       names <- as.numeric(unique(names))
  #       for(i in 1:length(names)){
  #         new[which(new == names[i])] <- input[[paste0("cluster_",names[i])]]
  #       }
  #     }else{
  #       names <- data[[input$assay.2]][[paste0('col_attrs/',input$group.by.2)]][]
  #       new <- names
  #       names <- unique(names)
  #       for(i in 1:length(names)){
  #         new[which(new == names[i])] <- input[[paste0(names[i],"_cells")]]
  #       }
  #     }
  #     new <- list(new)
  #     names(new) <- input$new_name
  #     if(!grepl("cell_type$",names(new))) names(new) <- paste0(names(new),"_cell_type")
  #     
  #     data[[input$assay.2]]$add.col.attribute(attributes = new, overwrite = TRUE)
  #     annot.trigger$trigger()
  #   }
  # })
  # 
  # observeEvent(input$add_to_tmp, {
  #   if(is.null(cells)){
  #     showNotification("Please select cells from either plot first", type = 'warning')
  #   }else if(is.null(input$custom_name) | input$custom_name == ""){
  #     showNotification("Please enter an annotation name", type = 'warning')
  #   }else{
  #     tmp_annotations[which(names(tmp_annotations)%in%cells[,1])] <<- input$custom_name
  #   }
  #   tmp.trigger$trigger()
  # })
  # 
  # observeEvent(input$save_custom_scheme, {
  #   if(all(tmp_annotations=='unlabled')){
  #     showNotification("No annotations have been added", type = 'warning')
  #   }else if(is.null(input$custom_scheme) | input$custom_scheme==""){
  #     showNotification("Please enter an name for the annotation scheme", type = 'warning')
  #   }else{
  #     new <- tmp_annotations
  #     names(new) <- NULL
  #     new <- list(tmp_annotations)
  #     names(new) <- input$custom_scheme
  #     if(!grepl("cell_type$",names(new))) names(new) <- paste0(names(new),"_cell_type")
  #     data[[input$assay.2]]$add.col.attribute(attributes = new, overwrite = TRUE)
  #     tmp_annotations <<- rep("unlabled",data[[input$assay.2]]$shape[2])
  #     annot.trigger$trigger()
  #   }
  # })
  # 
  # output$assay.2 <- renderUI({
  #   selectInput('assay.2', "Select Assay", choices = names(data), selected = ifelse(any(names(data)=="RNA"),yes = "RNA",no = names(data)[1]))
  # })
  # 
  # output$reduc.2 <- renderUI({
  #   req(input$assay.2)
  #   options <- c("tsne_rna_2d", "umap_rna_2d")
  #   meta.names <- data[[input$assay.2]][['col_attrs']]$names
  #   meta.names <- meta.names[sub("_.$","",meta.names)%in%options]
  #   meta.names <- unique(toupper(sub("_.*","",meta.names)))
  #   names(meta.names) <- meta.names
  #   selectInput(
  #     'reduc.2',
  #     "Select Reduction",
  #     choices = as.list(meta.names),
  #     selected = meta.names[1]
  #   )
  # })
  # 
  # cell_type_options_old <- NULL
  # output$group.by.2 <- renderUI({
  #   annot.trigger$depend()
  #   req(input$assay.2)
  #   
  #   cell_type_options <- names(data[[input$assay.2]]$col.attrs)[grep('cell_type|cell types|cell type|cell_types',names(data[[input$assay.2]]$col.attrs))]
  #   if(!is.null(cell_type_options_old) & (length(cell_type_options)!=length(cell_type_options_old))){
  #     sel <- cell_type_options[!cell_type_options%in%cell_type_options_old]
  #   }else{
  #     sel <- 'cluster'
  #   }
  #   cell_type_options_old <<- names(data[[input$assay.2]]$col.attrs)[grep('cell_type|cell types|cell type|cell_types',names(data[[input$assay.2]]$col.attrs))]
  #   selectInput(inputId = 'group.by.2', label = 'Choose Annotations', choices = c('cluster', cell_type_options), selected = sel, multiple = FALSE)
  # })
  # 
  # output$annots <- renderUI({
  #   annot.trigger$depend()
  #   req(input$group.by.2)
  # 
  #   if(input$group.by.2 == "cluster"){
  #     names <- as.numeric(unique(data[[input$assay.2]][[paste0("col_attrs/",tolower(input$assay.2),"_clusters")]][]))
  #     names <- names[order(names)]
  #     id <- "cluster"
  #     lapply(1:length(names), function(x){
  #       textInput(inputId = paste0(id,"_",names[x]), label = paste0(id," ",names[x]), value = paste0(id," ",names[x]))
  #     })
  #   }else{
  #     names <- unique(data[[input$assay.2]][[paste0('col_attrs/',input$group.by.2)]][])
  #     id <- "cells"
  #     lapply(1:length(names), function(x){
  #       if(!grepl("cells$",names[x])){
  #         textInput(inputId = paste0(names[x],"_cells"), label = paste0(names[x]," cells"), value = paste0(names[x]," cells"))
  #       }else{
  #         textInput(inputId = names[x], label = names[x], value = names[x])
  #       }
  #     })
  #   }
  # })
  # 
  # output$feature.select3 <- renderUI({
  #   req(input$assay.2)
  #   selectInput(inputId = 'feature.select3', label = 'Choose a Feature for Feature Plot', choices = data[[input$assay.2]][['row_attrs/features']][], selected = data[[input$assay.2]][['row_attrs/features']][1], multiple = FALSE)
  # })
  # 
  # meta.trigger <- makeReactiveTrigger()
  # 
  # order <- NULL
  # observeEvent(input$meta_group_checkbox, ignoreNULL = T,{
  #   new <- input$meta_group_checkbox
  #   #print(new)
  #   #print(order)
  #   if(length(new)>length(order) || is.null(order)){
  #     new <- new[!new%in%order]
  #     #print(new)
  #     order <<- c(order, new)
  #     #print(order)
  #   }else if(identical(new, character(0))){
  #     order <<- NULL
  #   }else{
  #     order <<- order[-which(!order%in%new)]
  #   }
  #   meta.trigger$trigger()
  # })
  # 
  # output$meta_group_checkbox <- renderUI({
  #   req(input$assay.3)
  #   checkboxGroupInput(inputId = 'meta_group_checkbox', label = 'Meta Data', choices = names(data[[input$assay.3]]$col.attrs))
  # })
  # output$example_meta_group_name <- renderText({
  #   meta.trigger$depend()
  #   if(length(order)>=2){
  #     paste(order,collapse=" + ")
  #   }else{
  #     NULL
  #   }
  # })
  # output$example_meta_group <- renderText({
  #   req(input$assay.3)
  #   meta.trigger$depend()
  #   if(length(order)>=2){
  #     ex <- NULL
  #     for(i in 1:length(order)){
  #       ex <- c(ex,data[[input$assay.3]][[paste0('col_attrs/',order[i])]][1])
  #     }
  #     paste(ex, collapse = "_")
  #   }else{
  #     NULL
  #   }
  # })
  # output$assay.3 <- renderUI({
  #   selectInput('assay.3', "Select Assay", choices = names(data), selected = ifelse(any(names(data)=="RNA"),yes = "RNA",no = names(data)[1]))
  # })
}

shinyApp(ui = ui, server = server)
