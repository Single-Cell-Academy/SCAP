# Server

server <- function(input, output, session){
  
  ###### set working dir and source required functions #####
  options(shiny.maxRequestSize=500*1024^2)
  setwd("/Users/jsjoyal/Desktop/SCAP/")
  volumes <- c(data = './test_data/')
  source("./R/SCAP_functions.R")
  
  ###### import data ######
  output$data_used <- renderText({
    if(!is.null(input$assay_1) & !is.na(input$directory[[1]][2])){
      return(paste0("Chosen data: ", input$directory[[1]][2]))
    }else{
      return(paste0("Please select your data..."))
    }
  })
  
  shinyDirChoose(input, "directory", roots = volumes, session = session, restrictions = system.file(package = "base"))
  
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
  
  ####################### MAIN TAB ########################
  #                                                       #
  ####################### MAIN TAB ########################
  
  # Golbal variables used by this tab
  grouping_options_old <- NULL
  
  #data_trigger <- makeReactiveTrigger()
  
  ########## Select Input for Assay ##############
  output$assay_1 <- renderUI({
    req(data_update())
    if(is.null(data)){
      return(NULL)
    }
    selectInput('assay_1', "Select Assay", choices = names(data), selected = ifelse(any(names(data)=="RNA"),yes = "RNA",no = names(data)[1]))
  })
  
  output$grouping_1 <- renderUI({
    req(data_update())
    req(input$assay_1)
    #data_trigger$depend()
    assay <- input$assay_1
    cols <- names(data[[assay]]$col.attrs)
    grouping_options <- cols[grep('_meta_data$',cols)]
    #if(!is.null(grouping_options_old) & (length(grouping_options)!= length(grouping_options_old))){
    #sel <- sub("_meta_data$", "", grouping_options[!grouping_options%in%grouping_options_old])
    if(any(grepl(paste0(tolower(assay),"_clusters"),sub("_meta_data","",grouping_options), fixed = TRUE))){
      sel <- paste0(tolower(assay),"_clusters")
    }else{
      sel <- sub("_meta_data","",grouping_options[1])
    }
    grouping_options_old <<- grouping_options
    
    selectInput(inputId = 'grouping_1', label = 'Group By', choices = sub("_meta_data$","",grouping_options), selected = sel, multiple = FALSE)
  })
  
  output$reduction_1 <- renderUI({
    req(data_update())
    req(input$assay_1)
    #data_trigger$depend()
    assay <- input$assay_1
    options <- names(data[[assay]]$col.attrs)
    l.assay <- tolower(assay)
    reductions <- unique(sub("_._reduction$","",options[which(grepl('_reduction$',options)&grepl(l.assay,options))]))
    selectInput(inputId = 'reduction_1', 'Choose Clustering Method', choices = as.list(reductions), selected = ifelse(test = any(reductions==paste0('tsne_',l.assay,'_2d')), yes = paste0('tsne_',l.assay,'_2d'), no = reductions[1]))
  })
  
  output$featureplot_1_feature.select<- renderUI({
    req(data_update())
    req(input$assay_1)
    #data_trigger$depend()
    assay <- input$assay_1
    selectInput(inputId = 'featureplot_1_feature.select', label = 'Select a Feature to Visualize on the Feature Plot (right)', choices = data[[assay]][['row_attrs/features']][], selected = data[[assay]][['row_attrs/features']][1], multiple = FALSE)
  })
  output$dotplot_1_feature_select<- renderUI({
    req(data_update())
    req(input$assay_1)
    #data_trigger$depend()
    assay <- input$assay_1
    selectInput(inputId = 'dotplot_1_feature_select', label = 'Choose Features to Visualize on Dot Plot (bottom)', choices = data[[assay]][['row_attrs/features']][], selected = data[[assay]][['row_attrs/features']][3], multiple = TRUE)
  })
  
  ############### Dim Plot ################
  output$dimplot_1 <- renderPlotly({
    req(data_update())
    req(input$grouping_1)
    dimPlotlyOutput(assay.in = input$assay_1, reduc.in = input$reduction_1, group.by = input$grouping_1, data = data)
  })
  
  ########## Feature Plot #####################
  output$featureplot_1 <- renderPlotly({
    req(data_update())
    req(input$featureplot_1_feature.select)
    featurePlotlyOutput(assay.in = input$assay_1, reduc.in = input$reduction_1, group.by = input$grouping_1, feature.in = input$featureplot_1_feature.select, data = data)
  })
  
  ######### Dot Plot #########
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
  
  
  ####################### ANNOTATION TAB ########################
  #                                                             #
  ####################### ANNOTATION TAB ########################
  
  # Global variables used by this tab
  tmp_annotations <- NULL
  meta_options_old <- NULL
  cells <- NULL
  
  # Reactive triggers used by this tab
  tmp.trigger <- makeReactiveTrigger()
  annot.trigger <- makeReactiveTrigger()
  
  output$assay_2 <- renderUI({
    req(data_update())
    #data_trigger$depend()
    selectInput('assay_2', "Select Assay", choices = names(data), selected = ifelse(any(names(data)=="RNA"),yes = "RNA",no = names(data)[1]))
  })
  
  output$grouping_2 <- renderUI({
    req(data_update())
    annot.trigger$depend()
    req(input$assay_2)
    #data_trigger$depend()
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
  
  output$reduction_2 <- renderUI({
    req(data_update())
    req(input$assay_2)
    #data_trigger$depend()
    assay <- input$assay_2
    options <- names(data[[assay]]$col.attrs)[grep("_2d",names(data[[assay]]$col.attrs))]
    l.assay <- tolower(assay)
    reductions <- unique(sub("_._reduction$","",options[which(grepl('_reduction$',options)&grepl(l.assay,options))]))
    selectInput(inputId = 'reduction_2', 'Choose Clustering Method', choices = as.list(reductions), selected = ifelse(test = any(reductions==paste0('tsne_',l.assay,'_2d')), yes = paste0('tsne_',l.assay,'_2d'), no = reductions[1]))
  })
  
  output$featureplot_2_feature_select<- renderUI({
    req(data_update())
    req(input$assay_2)
    #data_trigger$depend()
    assay <- input$assay_2
    selectInput(inputId = 'featureplot_2_feature_select', label = 'Select a Feature to Visualize on the Feature Plot (right)', choices = data[[assay]][['row_attrs/features']][], selected = data[[assay]][['row_attrs/features']][1], multiple = FALSE)
  })
  
  output$selected_cells <- renderTable({
    req(data_update())
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
  
  output$dimplot_2 <- renderPlotly({
    req(data_update())
    req(input$grouping_2)
    dimPlotlyOutput(assay.in = input$assay_2, reduc.in = input$reduction_2, group.by = input$grouping_2, data = data)
  })
  
  output$featureplot_2 <- renderPlotly({
    req(data_update())
    req(input$featureplot_2_feature_select)
    featurePlotlyOutput(assay.in = input$assay_2, reduc.in = input$reduction_2, group.by = input$grouping_2, feature.in = input$featureplot_2_feature_select, data = data)
  })
  
  observeEvent(input$find.markers, ignoreNULL = FALSE, ignoreInit = FALSE, handlerExpr = {
    output$markers <- renderTable({
      req(input$grouping_2)
      #data_trigger$depend()
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
  
  
  
  #---------------------------------------------------#
  
  # output$annot.select.cluster <- renderPlotly({
  #   tmp.trigger$depend()
  #   annot.trigger$depend()
  #   req(input$reduc.2)
  #   req(input$group.by.2)
  #   req(input$assay_2)
  #   
  #   reduc <- tolower(paste0(input$reduc.2,"_rna_2d_",1:2))
  #   
  #   if(input$group.by.2 == "cluster"){
  #     plot.data <- data[[input$assay_2]]$get.attribute.df(attributes=c(reduc,tolower(paste0(input$assay_2,"_","clusters")),'percent.mt'))
  #   }
  #   else{
  #     plot.data <- data[[input$assay_2]]$get.attribute.df(attributes=c(reduc,input$group.by.2,'percent.mt'))
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
  #   ) %>% layout(title = ifelse(test = input$annot_panel == 'cluster', yes = paste0(input$assay_2, " data coloured by ", input$group.by.2), no = paste0(input$assay_2, " data")),xaxis = ax.x, yaxis = ax.y, dragmode = "select")
  # 
  #   p
  # })
  # output$annot.select.feature <- renderPlotly({
  #   tmp.trigger$depend()
  #   annot.trigger$depend()
  #   req(input$feature.select3)
  #   req(input$reduc.2)
  #   req(input$assay_2)
  #   data.features <- as.data.frame(data[[input$assay_2]][['matrix']][,which(data[[input$assay_2]]$row.attrs$features[]%in%input$feature.select3)])
  #   if(nrow(data.features)==0|ncol(data.features)==0) return(NULL)
  #   if(input$group.by.2 == 'cluster'){
  #     data.annot <- data[[input$assay_2]]$get.attribute.df(attributes=tolower(paste0(input$assay_2,"_clusters")))
  #   }else{
  #     data.annot <- data[[input$assay_2]]$get.attribute.df(attributes=input$group.by.2)
  #   }
  #   rownames(data.features) <- data[[input$assay_2]]$col.attrs$CellID[]
  #   colnames(data.features) <- input$feature.select3
  # 
  #   
  #   dims <- tolower(paste0(input$reduc.2,"_",input$assay_2,"_2d_",1:2))
  # 
  #   plot.data <- data[[input$assay_2]]$get.attribute.df(attributes=dims)
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
  #     req(input$assay_2)
  #     req(input$group.by.2)
  #     #if(is.null(input$find.markers)) NULL
  #     m <- t(data[[input$assay_2]]$matrix[,])
  #     rownames(m) <- data[[input$assay_2]]$row.attrs$features[]
  #     colnames(m) <- data[[input$assay_2]]$col.attrs$CellID[]
  #     if(input$group.by.2=='cluster'){
  #       y <- data[[input$assay_2]][[paste0("col_attrs/",tolower(input$assay_2),"_clusters")]][]
  #     }else{
  #       y <- data[[input$assay_2]][[paste0('col_attrs/',input$group.by.2)]][]
  #     }
  #     cols <- unlist(lapply(cells[,1], function(x){
  #       grep(x,data[[input$assay_2]]$col.attrs$CellID[])
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
  #   req(input$assay_2)
  #   req(input$group.by.2)
  #   if(is.null(input$new_name) || input$new_name == "" || input$new_name == " "){
  #     showNotification("Please name the annotation scheme", type = 'warning')
  #   }else{
  #     if(input$group.by.2 == "cluster"){
  #       names <- data[[input$assay_2]][[paste0("col_attrs/",tolower(input$assay_2),"_clusters")]][]
  #       new <- names
  #       names <- as.numeric(unique(names))
  #       for(i in 1:length(names)){
  #         new[which(new == names[i])] <- input[[paste0("cluster_",names[i])]]
  #       }
  #     }else{
  #       names <- data[[input$assay_2]][[paste0('col_attrs/',input$group.by.2)]][]
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
  #     data[[input$assay_2]]$add.col.attribute(attributes = new, overwrite = TRUE)
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
  #     data[[input$assay_2]]$add.col.attribute(attributes = new, overwrite = TRUE)
  #     tmp_annotations <<- rep("unlabled",data[[input$assay_2]]$shape[2])
  #     annot.trigger$trigger()
  #   }
  # })
  # 
  # output$assay_2 <- renderUI({
  #   selectInput('assay_2', "Select Assay", choices = names(data), selected = ifelse(any(names(data)=="RNA"),yes = "RNA",no = names(data)[1]))
  # })
  # 
  # output$reduc.2 <- renderUI({
  #   req(input$assay_2)
  #   options <- c("tsne_rna_2d", "umap_rna_2d")
  #   meta.names <- data[[input$assay_2]][['col_attrs']]$names
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
  #   req(input$assay_2)
  #   
  #   cell_type_options <- names(data[[input$assay_2]]$col.attrs)[grep('cell_type|cell types|cell type|cell_types',names(data[[input$assay_2]]$col.attrs))]
  #   if(!is.null(cell_type_options_old) & (length(cell_type_options)!=length(cell_type_options_old))){
  #     sel <- cell_type_options[!cell_type_options%in%cell_type_options_old]
  #   }else{
  #     sel <- 'cluster'
  #   }
  #   cell_type_options_old <<- names(data[[input$assay_2]]$col.attrs)[grep('cell_type|cell types|cell type|cell_types',names(data[[input$assay_2]]$col.attrs))]
  #   selectInput(inputId = 'group.by.2', label = 'Choose Annotations', choices = c('cluster', cell_type_options), selected = sel, multiple = FALSE)
  # })
  # 
  # output$annots <- renderUI({
  #   annot.trigger$depend()
  #   req(input$group.by.2)
  # 
  #   if(input$group.by.2 == "cluster"){
  #     names <- as.numeric(unique(data[[input$assay_2]][[paste0("col_attrs/",tolower(input$assay_2),"_clusters")]][]))
  #     names <- names[order(names)]
  #     id <- "cluster"
  #     lapply(1:length(names), function(x){
  #       textInput(inputId = paste0(id,"_",names[x]), label = paste0(id," ",names[x]), value = paste0(id," ",names[x]))
  #     })
  #   }else{
  #     names <- unique(data[[input$assay_2]][[paste0('col_attrs/',input$group.by.2)]][])
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
  #   req(input$assay_2)
  #   selectInput(inputId = 'feature.select3', label = 'Choose a Feature for Feature Plot', choices = data[[input$assay_2]][['row_attrs/features']][], selected = data[[input$assay_2]][['row_attrs/features']][1], multiple = FALSE)
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