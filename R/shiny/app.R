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
if(!require('loomR'))
  BiocManager::install('loomR')
library(loomR)



#---------------------------------------------------------#
#                   User Interface                        #
#_________________________________________________________#

ui <- navbarPage(fluid = TRUE,
  theme = shinytheme('yeti'),
  title = "Single Cell Analysis Portal",
  tabPanel(
    "Main",
    sidebarPanel(
      "",
      uiOutput('assay'),
      uiOutput("reduc"),
      uiOutput('feature.select'),
      uiOutput('feature.select2'),
      uiOutput('group.by')
    ),
    mainPanel(
      "",
      fluidRow(style='padding:100px;',
        column(
          width = 6,
          plotlyOutput('cluster', height = '500px', width = '500px')%>% withSpinner(color="black")
        ),
        column(
          width = 6,
          plotlyOutput('feature', height = '500px', width = '500px')%>% withSpinner(color="black")
        ),
        div(style = "height:10px;")
      ),
      fluidRow(
        column(width = 12, offset = 1,
               plotOutput('dp', height = '500px', width = '1000px')%>% withSpinner(color="black")
        )
      )
    )
  ),
  tabPanel(
    "Cell Annotation",
    sidebarPanel(
      uiOutput('assay.2'),
      uiOutput('reduc.2'),
      uiOutput('feature.select3'),
      uiOutput('group.by.2'),
      uiOutput('annots'),
      h4("Selected Cells"),
      tableOutput('cells'),
      tags$head(tags$style("#cells{overflow-y:scroll; max-height: 200px; background: ghostwhite;}"))
    ),
    mainPanel(
      fluidRow(style='padding:100px;',
        column(width = 6,
               plotlyOutput('annot.select.cluster', height = '500px', width = '500px')%>% withSpinner(color="black")
               ),
        column(width = 6,
               plotlyOutput('annot.select.feature', height = '500px', width = '500px')%>% withSpinner(color="black")
               )
      ),
      fluidRow(
        column(width = 12,style='padding:10px;',
               actionButton('find.markers', "Find Markers for Selected Cells"),
               tableOutput('markers') %>% withSpinner(color="black")
               )
      )
    )
  )
)
#jqui_resizable(jqui_draggable(
#%>% withSpinner(color="#0dc5c1")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


#---------------------------------------------------------#
#                       SERVER                            #
#_________________________________________________________#

server <- function(input, output){

  ###### set working dir and source required functions #####
  options(shiny.maxRequestSize=500*1024^2)
  setwd("/Users/jsjoyal/Desktop/SCAP/")
  source("./R/SCAP_functions.R")
  
  ###### import data ######
  files <- list.files("./test_data/test_project/")
  assays <- assays <- sub(".loom","",files)
  data <- list()
  for(i in 1:length(assays)){
    data[[i]] <- connect(paste0("./test_data/test_project/",files[i]), mode = "r+")
  }
  names(data) <- assays
  
  ####################### MAIN TAB ########################
  #                                                       #
  ####################### MAIN TAB ########################
  
  ######### Dot Plot #########
  output$dp <- renderPlot({
    req(input$assay)
    req(input$feature.select)
    if(input$group.by == 'cluster'){
      group.by <- tolower(paste0(input$assay,"_clusters"))
    }else{
      group.by <- "cell_type"
    }
    
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
    
    #p <- ggplotly(dotPlot(data = data[[input$assay]], assay = input$assay, features = input$feature.select, group.by = group.by))
    p <- dotPlot(data = data[[input$assay]], assay = input$assay, features = input$feature.select, group.by = group.by) + guides(color = guide_colourbar(order = 1, title = "Average Expression"), size = guide_legend(order = 2, title = "Percent Expressed"))
    if(identical(class(p),"shiny.tag")) return(NULL)
    #p$x$data[[1]]$hovertemplate <- paste("<b>pct.exp: %{~pct.exp}","<extra></extra>")
    #p %>% layout(xaxis = ax.x, yaxis = ax.y,plot_bgcolor="white")
    p + theme_few()
  })
 
  ############### Dim Plot ################
  output$cluster <- renderPlotly({
    req(input$assay)
    req(input$reduc)
    req(input$group.by)
    
    n <- if(grepl("2D", input$reduc)){
      2
    }else{
      3
    }
    
    reduc <- tolower(paste0(sub(" ",paste0("_",input$assay,"_"),input$reduc),"_",1:n))
    
    if(input$group.by == "cluster"){
      plot.data <- data[[input$assay]]$get.attribute.df(attributes=c(reduc,tolower(paste0(input$assay,"_","clusters")),'percent.mt'))
    }
    else{
      plot.data <- data[[input$assay]]$get.attribute.df(attributes=c(reduc,'cell_type','percent.mt'))
    }
    key <- reduc_key(key = toupper(sub(" ..","",input$reduc)))
    
    ax.x <- list(
      title = paste0(key,"_1"),
      zeroline = FALSE,
      showline = TRUE,
      showticklabels = FALSE,
      showgrid = FALSE,
      mirror=TRUE,
      ticks='none'
    )
    ax.y <- list(
      title = paste0(key,"_2"),
      zeroline = FALSE,
      showline = TRUE,
      showticklabels = FALSE,
      showgrid = FALSE,
      mirror=TRUE,
      ticks='none'
    )
    
    ax.z <- list(
      title = paste0(key,"_3"),
      zeroline = FALSE,
      showline = TRUE,
      showticklabels = FALSE,
      showgrid = FALSE,
      mirror=TRUE,
      ticks='none'
    )
    
    #print(head(plot.data))
    if(n == 2){
      p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2], 
              type = 'scatter', mode = 'markers', 
              color = plot.data[,3], text =  ~paste0(
                key,"_1: ", format(plot.data[,1],digits=3),"\n",
                "</br>",key,"_2: ", format(plot.data[,2],digits=3), "\n",
                "</br>percent.mt: ", format(plot.data[,4],digits=3), "%"), 
              hovertemplate = paste0('<b>%{text}</b>')
              ) %>% layout(title = paste0(input$assay, " data coloured by ", input$group.by) ,xaxis = ax.x, yaxis = ax.y)
    }else if(n==3){
      p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2], z = plot.data[,3],
                   type = 'scatter3d', mode = 'markers', 
                   color = plot.data[,4], text =  ~paste0(
                     key,"_1: ", format(plot.data[,1],digits=3),"\n",
                     "</br>",key,"_2: ", format(plot.data[,2],digits=3), "\n",
                     "</br>",key,"_3: ", format(plot.data[,3],digits=3), "\n",
                     "</br>percent.mt: ", format(plot.data[,5],digits=3), "%"), 
                   hovertemplate = paste0('<b>%{text}</b>')
                   ) %>% layout(title = paste0(input$assay, " data coloured by ", input$group.by) ,scene = list(xaxis = ax.x, yaxis = ax.y, zaxis = ax.z))
    }else NULL
    p
  })
  
  ########## Feature Plot #####################
  output$feature <- renderPlotly({
    req(input$assay)
    req(input$feature.select2)
    req(input$reduc)
    data.features <- as.data.frame(data[[input$assay]][['matrix']][,which(data[[input$assay]]$row.attrs$features[]%in%input$feature.select2)])
    if(nrow(data.features)==0|ncol(data.features)==0) return(NULL)
    if(input$group.by == 'cluster'){
      data.annot <- data[[input$assay]]$get.attribute.df(attributes=tolower(paste0(input$assay,"_clusters")))
    }else{
      data.annot <- data[[input$assay]]$get.attribute.df(attributes='cell_type')
    }
    rownames(data.features) <- data[[input$assay]]$col.attrs$CellID[]
    colnames(data.features) <- input$feature.select2
    
    n <- if(grepl("2D", input$reduc)){
      2
    }else{
      3
    }
    
    dims <- tolower(paste0(sub(" ",paste0("_",input$assay,"_"),input$reduc),"_",1:n))
    
    plot.data <- data[[input$assay]]$get.attribute.df(attributes=dims)
    plot.data <- cbind(plot.data,data.annot)
    plot.data <- cbind(plot.data,data.features)
    
    #plot.data[,4] <- scale(plot.data[,4])
    
    key <- reduc_key(key = toupper(sub(" ..","",input$reduc)))
    
    ax.x <- list(
      title = paste0(key,"_1"),
      zeroline = FALSE,
      showline = TRUE,
      showticklabels = FALSE,
      showgrid = FALSE,
      mirror=TRUE,
      ticks='none'
    )
    ax.y <- list(
      title = paste0(key,"_2"),
      zeroline = FALSE,
      showline = TRUE,
      showticklabels = FALSE,
      showgrid = FALSE,
      mirror=TRUE,
      ticks='none'
    )
    ax.z <- list(
      title = paste0(key,"_3"),
      zeroline = FALSE,
      showline = TRUE,
      showticklabels = FALSE,
      showgrid = FALSE,
      mirror=TRUE,
      ticks='none'
    )
    
    if(n == 2){
      plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2],
              type = 'scatter', mode = 'markers', 
              color = plot.data[,4],text =  ~paste0(
                key,"_1: ", format(plot.data[,1],digits=3),"\n",
                "</br>",key,"_2: ", format(plot.data[,2],digits=3), "\n",
                "</br>",input$group.by,": ", plot.data[,3]), 
              hovertemplate = paste0('<b>%{text}</b>',
                                     '<extra></extra>')
      ) %>% layout(title = input$feature.select2 ,xaxis = ax.x, yaxis = ax.y)
    }else if(n==3){
      plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2], z = plot.data[,3],
              type = 'scatter3d', mode = 'markers', 
              color = plot.data[,5],text =  ~paste0(
                key,"_1: ", format(plot.data[,1],digits=3),"\n",
                "</br>",key,"_2: ", format(plot.data[,2],digits=3), "\n",
                "</br>",key,"_3: ", format(plot.data[,3],digits=3), "\n",
                "</br>",input$group.by,": ", plot.data[,4]), 
              hovertemplate = paste0('<b>%{text}</b>',
                                     '<extra></extra>')
      ) %>% layout(title = input$feature.select2 , scene = list(xaxis = ax.x, yaxis = ax.y, zaxis = ax.z))
    }else NULL
  })
  
  ########## Select Input for Assay ##############
  output$assay <- renderUI({
    selectInput('assay', "Select Assay", choices = names(data), selected = ifelse(any(names(data)=="RNA"),yes = "RNA",no = names(data)[1]))
  })
  
  ######### Select input for reduction ##########
  output$reduc <- renderUI({
    req(input$assay)
    options <- c("pca", "tsne", "umap", "diff")
    meta.names <- data[[input$assay]][['col_attrs']]$names
    meta.names <- unique(sub(paste0(tolower(input$assay),"_"),"",sub("_.$","",meta.names[grep(paste(options,collapse = "|"),meta.names)])))
    meta.names <- toupper(sub("_"," ",meta.names))
    names(meta.names) <- meta.names
    selectInput(
      'reduc',
      "Select Reduction",
      choices = as.list(meta.names),
      selected = ifelse(any(meta.names == "TSNE 2D"), yes = "TSNE 2D", no = meta.names[1])
    )
  })
 
  ########## Select features for dotplot ###########
  output$feature.select <- renderUI({
    req(input$assay)
    selectInput(inputId = 'feature.select', label = 'Choose Features for Dot Plot', choices = data[[input$assay]][['row_attrs/features']][], selected = data[[input$assay]][['row_attrs/features']][1:3], multiple = TRUE)
  })

  ########## Select features for feature plot ###########
  output$feature.select2 <- renderUI({
    req(input$assay)
    selectInput(inputId = 'feature.select2', label = 'Choose a Feature for Feature Plot', choices = data[[input$assay]][['row_attrs/features']][], selected = data[[input$assay]][['row_attrs/features']][1], multiple = FALSE)
  })
  
  ########## Select group.by ###########
  output$group.by <- renderUI({
    selectInput(inputId = 'group.by', label = 'Choose Annotations', choices = c('cluster', 'cell type'), selected = 'cluster', multiple = FALSE)
  })
 
  ####################### ANNOTATION TAB ########################
  #                                                             #
  ####################### ANNOTATION TAB ########################
  
  output$annot.select.cluster <- renderPlotly({
      req(input$reduc.2)
      req(input$group.by.2)
      req(input$assay.2)
      
      reduc <- tolower(paste0(input$reduc.2,"_rna_2d_",1:2))
      
      if(input$group.by.2 == "cluster"){
        plot.data <- data[[input$assay.2]]$get.attribute.df(attributes=c(reduc,tolower(paste0(input$assay.2,"_","clusters")),'percent.mt'))
      }
      else{
        plot.data <- data[[input$assay.2]]$get.attribute.df(attributes=c(reduc,'cell_type','percent.mt'))
      }
      key <- reduc_key(key = toupper(input$reduc.2))
      
      ax.x <- list(
        title = paste0(key,"_1"),
        zeroline = FALSE,
        showline = TRUE,
        showticklabels = FALSE,
        showgrid = FALSE,
        mirror=TRUE,
        ticks='none'
      )
      ax.y <- list(
        title = paste0(key,"_2"),
        zeroline = FALSE,
        showline = TRUE,
        showticklabels = FALSE,
        showgrid = FALSE,
        mirror=TRUE,
        ticks='none'
      )
      
      ax.z <- list(
        title = paste0(key,"_3"),
        zeroline = FALSE,
        showline = TRUE,
        showticklabels = FALSE,
        showgrid = FALSE,
        mirror=TRUE,
        ticks='none'
      )
      
      #print(head(plot.data))
      
      p <- plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2], key = ~rownames(plot.data),
                   type = 'scatter', mode = 'markers', 
                   color = plot.data[,3], text =  ~paste0(
                     key,"_1: ", format(plot.data[,1],digits=3),"\n",
                     "</br>",key,"_2: ", format(plot.data[,2],digits=3), "\n",
                     "</br>percent.mt: ", format(plot.data[,4],digits=3), "%"), 
                   hovertemplate = paste0('<b>%{text}</b>')
      ) %>% layout(title = paste0(input$assay.2, " data coloured by ", input$group.by.2) ,xaxis = ax.x, yaxis = ax.y, dragmode = "select")
      
      p
  })
  output$annot.select.feature <- renderPlotly({
    req(input$feature.select3)
    req(input$reduc.2)
    req(input$assay.2)
    data.features <- as.data.frame(data[[input$assay.2]][['matrix']][,which(data[[input$assay.2]]$row.attrs$features[]%in%input$feature.select3)])
    if(nrow(data.features)==0|ncol(data.features)==0) return(NULL)
    if(input$group.by.2 == 'cluster'){
      data.annot <- data[[input$assay.2]]$get.attribute.df(attributes=tolower(paste0(input$assay.2,"_clusters")))
    }else{
      data.annot <- data[[input$assay.2]]$get.attribute.df(attributes='cell_type')
    }
    rownames(data.features) <- data[[input$assay.2]]$col.attrs$CellID[]
    colnames(data.features) <- input$feature.select3

    
    dims <- tolower(paste0(input$reduc.2,"_",input$assay.2,"_2d_",1:2))

    plot.data <- data[[input$assay.2]]$get.attribute.df(attributes=dims)
    plot.data <- cbind(plot.data,data.annot)
    plot.data <- cbind(plot.data,data.features)

    #plot.data[,4] <- scale(plot.data[,4])

    key <- reduc_key(key = toupper(input$reduc.2))

    ax.x <- list(
      title = paste0(key,"_1"),
      zeroline = FALSE,
      showline = TRUE,
      showticklabels = FALSE,
      showgrid = FALSE,
      mirror=TRUE,
      ticks='none'
    )
    ax.y <- list(
      title = paste0(key,"_2"),
      zeroline = FALSE,
      showline = TRUE,
      showticklabels = FALSE,
      showgrid = FALSE,
      mirror=TRUE,
      ticks='none'
    )
    ax.z <- list(
      title = paste0(key,"_3"),
      zeroline = FALSE,
      showline = TRUE,
      showticklabels = FALSE,
      showgrid = FALSE,
      mirror=TRUE,
      ticks='none'
    )

    
    plot_ly(data = plot.data, x = plot.data[,1], y = plot.data[,2],
            type = 'scatter', mode = 'markers', key = ~rownames(plot.data),
            color = plot.data[,4],text =  ~paste0(
              key,"_1: ", format(plot.data[,1],digits=3),"\n",
              "</br>",key,"_2: ", format(plot.data[,2],digits=3), "\n",
              "</br>",input$group.by.2,": ", plot.data[,3]),
            hovertemplate = paste0('<b>%{text}</b>',
                                   '<extra></extra>')
    ) %>% layout(title = input$feature.select3 ,xaxis = ax.x, yaxis = ax.y, dragmode = "select")
  })
  
  cells <- NULL
  output$cells <- renderTable({
    d <- as.data.frame(event_data("plotly_selected")$key, stringsAsFactors = FALSE)
    if (is.null(d) | ncol(d) == 0 | nrow(d) == 0){
      cells <<- NULL
      "Click and drag to select cells for marker identification (i.e., select/lasso) (double-click to clear)"
      }else{
        colnames(d) <- ""
        #print(d)
        cells <<- d
        d
      }
  })
  
  observeEvent(input$find.markers, ignoreNULL = FALSE, ignoreInit = FALSE, handlerExpr = {
    output$markers <- renderTable({
      req(input$assay.2)
      #if(is.null(input$find.markers)) NULL
      m <- t(data[[input$assay.2]]$matrix[,])
      rownames(m) <- data[[input$assay.2]]$row.attrs$features[]
      colnames(m) <- data[[input$assay.2]]$col.attrs$CellID[]
      if(input$group.by.2=='cluster'){
        y <- data[[input$assay.2]][[paste0("col_attrs/",tolower(input$assay.2),"_clusters")]][]
      }else{
        y <- data[[input$assay.2]]$col.attrs$cell_type[]
      }
      cols <- unlist(lapply(cells[,1], function(x){
        grep(x,data[[input$assay.2]]$col.attrs$CellID[])
      }))
      #print(cells)
      #print(cols)
      #print(head(cols))
      if(!is.null(cells))
        y[cols] <- 'Selected'
      #print(head(y))
      t <- as.data.frame(top_markers(wilcoxauc(X = m, y = y), n = 10))[1:10,]
      if(!is.null(cells)){
        # t1 <- t[,-which(colnames(t)=="Selected")]
        # t <- cbind(t1,t[,which(colnames(t)=="Selected")])
        # colnames(t)[ncol(t)] <- "Selected"
        t <- t[,which(colnames(t)=="Selected")]
      }
      t
    }, hover = TRUE, width = '100%')
  })
  
  output$assay.2 <- renderUI({
    selectInput('assay.2', "Select Assay", choices = names(data), selected = ifelse(any(names(data)=="RNA"),yes = "RNA",no = names(data)[1]))
  })

  output$reduc.2 <- renderUI({
    req(input$assay.2)
    options <- c("tsne_rna_2d", "umap_rna_2d")
    meta.names <- data[[input$assay.2]][['col_attrs']]$names
    meta.names <- meta.names[sub("_.$","",meta.names)%in%options]
    meta.names <- unique(toupper(sub("_.*","",meta.names)))
    names(meta.names) <- meta.names
    selectInput(
      'reduc.2',
      "Select Reduction",
      choices = as.list(meta.names),
      selected = meta.names[1]
    )
  })

  output$group.by.2 <- renderUI({
    selectInput(inputId = 'group.by.2', label = 'Choose Annotations', choices = c('cluster', 'cell type'), selected = 'cluster', multiple = FALSE)
  })

  output$annots <- renderUI({
    req(input$group.by.2)

    if(input$group.by.2 == "cluster"){
      names <- as.numeric(unique(data[[input$assay.2]][[paste0("col_attrs/",tolower(input$assay.2),"_clusters")]][]))
      names <- names[order(names)]
      id <- "cluster"
      lapply(1:length(names), function(x){
        textInput(inputId = paste0(id,"_",names[x]), label = paste0(id," ",names[x]), value = paste0(id," ",names[x]))
      })
    }else{
      names <- unique(data[[input$assay.2]][['col_attrs/cell_type']][])
      id <- "cells"
      lapply(1:length(names), function(x){
        textInput(inputId = paste0(names[x],"_cells"), label = paste0(names[x]," cells"), value = paste0(names[x]," cells"))
      })
    }
  })

  output$feature.select3 <- renderUI({
    req(input$assay.2)
    selectInput(inputId = 'feature.select3', label = 'Choose a Feature for Feature Plot', choices = data[[input$assay.2]][['row_attrs/features']][], selected = data[[input$assay.2]][['row_attrs/features']][1], multiple = FALSE)
  })
  
}

shinyApp(ui = ui, server = server)
