library("cowplot")
library("dplyr")
library("ggplot2")
library("ggthemes")
library("gtools")
library("Matrix")
library("MODIS")
library("plotly")
library("rjson")
library("shiny")
library("shinyFiles")
library("data.table")
library("scibet")
library("readr")
library("reactable")
library("reticulate")
library("shinyjs")
library("presto")
library("bbplot")

reticulate::use_virtualenv("../renv/python/virtualenvs/renv-python-3.8.5/")

#### Variables that persist across sessions
## Read in table with datasets available for SciBet
datasets_scibet <- fread("../meta/SciBet_reference_list.tsv")
## Source functions
source("SCAP_functions.R")
source_python("../Python/rank_genes_groups_df.py")

anndata <- import('anndata')
scanpy <- import('scanpy')

init <- 0 # flag for autosave

server <- function(input, output, session){
  session$onSessionEnded(stopApp)
  
  options(shiny.maxRequestSize=500*1024^2)

  rvalues <- reactiveValues(tmp_annotations = NULL, cells = NULL, order = NULL, features = NULL, obs = NULL, obs_cat = NULL, reductions = NULL, cell_ids = NULL, h5ad = NULL, path_to_data = NULL,
                            raw_dtype = NULL)
  rvalues_mod <- reactiveValues(tmp_annotations = NULL, cells = NULL, order = NULL, features = NULL, obs = NULL, obs_cat = NULL, reductions = NULL, cell_ids = NULL, h5ad = NULL, path_to_data = NULL,
                                raw_dtype = NULL)

  ## Determine folders for ShinyDir button
  volumes <- c("FTP" = "/ftp", Home = fs::path_home())
  
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
  shinyFileChoose(input, "h5ad_in", roots = volumes, session = session)
  
  # connect chosen .h5ad file
  observeEvent(input$h5ad_in, {
    path <- parseFilePaths(selection = input$h5ad_in, roots = volumes)$datapath
    if(is.integer(path[1]) || identical(path, character(0)) || identical(path, character(0))) return(NULL)
    h5ad_files <- path#paste0(path,"/",list.files(path))
    assays <- sub(".h5ad","",sub(paste0(".*/"),"",h5ad_files))
    data <- list()
    ## Iterate over all assays and connect to h5ad objects
    for(i in 1:length(assays)){
      data[[i]] <- tryCatch({
        anndata$read(h5ad_files[i])
      },
      error = function(e){
        showModal(modalDialog(p(paste0("An error occured trying to connect to ", h5ad_files[i])), title = "Error connecting to h5ad file."), session = getDefaultReactiveDomain())
        return(NULL)
      })
    }
    if(is.null(data)) return(NULL)
    if(length(data) != length(assays)) return(NULL)
    if(length(unlist(lapply(data, function(x){x}))) != length(assays)) return(NULL)
    names(data) <- assays
    ## Check if RAW Anndata object is present or not. If not present, use the main object
    if(is.null(data[[1]]$raw)){ 
      rvalues$features <- rownames(data[[1]]$var)
    }else{
      test_gene_name <- rownames(data[[1]]$var)[1]
      if(test_gene_name %in% rownames(data[[1]]$raw$var)){ # check if rownames are numbers or gene names
        rvalues$features <- rownames(data[[1]]$raw$var)
      }else if("features" %in% colnames(data[[1]]$raw$var)){ ## Check if there is a column named features in raw
        rvalues$features <- data[[1]]$raw$var$features
      }else if(test_gene_name %in% data[[1]]$raw$var[,1]){ # otherwise, check if the first column contains rownames
        rvalues$features <- data[[1]]$raw$var[,1]
      }
    }
    rvalues$obs <- data[[1]]$obs_keys()
    ## Determine type of annotation and create a layer to annotate for easy usage later on
    rvalues$obs_cat <- check_if_obs_cat(obs_df = data[[1]]$obs) ## Function to check if an observation is categorical or numeric
    reductions <- data[[1]]$obsm$as_dict()
    if(length(reductions) == 0){
        showModal(modalDialog(p(paste0(h5ad_files[i], " has no dimensional reductions.")), title = "Error connecting to h5ad file."), session = getDefaultReactiveDomain())
        return(NULL)
    }
    reduction_keys <- data[[1]]$obsm_keys()
    r_names <- rownames(data[[1]]$obs)
    for(i in 1:length(reductions)){
      reductions[[i]] <- as.data.frame(reductions[[i]])
      colnames(reductions[[i]]) <- paste0(reduction_keys[i], "_", 1:ncol(reductions[[i]]))
      rownames(reductions[[i]]) <- r_names
    }
    names(reductions) <- reduction_keys
    rvalues$reductions <- reductions
    rvalues$cell_ids <- rownames(data[[1]]$obs)
    rvalues$h5ad <- data
    rvalues$path_to_data <- h5ad_files
      
    ## unload modality rvalues
    for(i in names(rvalues_mod)){
      rvalues_mod[[i]] <- NULL
    }
      
    ## Determine what data is likely stored in .raw
    if(is.null(data[[1]]$raw)){ ## Check if raw exists
      rvalues$raw_dtype <- "NULL"
    }else if(sum(rvalues$h5ad[[1]]$raw$X[1,]) %% 1 == 0){ ## Check whether raw contains un-normalized data or normalized data
      rvalues$raw_dtype <- "counts"
    }else{ ## Only if the other two conditions fail, use raw values to calculate differential expression
      rvalues$raw_dtype <- "normalized"
    }
      
    init <<- 0
  })

  # observe({ # auto save h5ad file(s)
  #   req(rvalues$h5ad)
  #   invalidateLater(120000) # 2 min
  #   if(init>0){
  #     #tryCatch(
  #     #  { 
  #         cat(file = stderr(), paste0(rvalues$path_to_data, "\n"))
  #         showNotification("Saving...", duration = NULL, id = 'auto_save')
  #         for(i in 1:length(rvalues$path_to_data)){
  #           rvalues$h5ad[[i]]$write(filename = rvalues$path_to_data[i])
  #         }
  #         removeNotification(id = 'auto_save')
  #      # },
  #      # error = function(e)
  #      # {
  #         #cat(file = stderr(), unlist(e))
  #     #    showModal(modalDialog(p(paste0("An error occured trying to write to ", rvalues$path_to_data[i], ": ", unlist(e))), title = "Error writing to h5ad file."), session = getDefaultReactiveDomain())
  #     #  }
  #    # )
  #   }
  #   init <<- init + 1
  # })

  source(file.path("server", "main.server.R"),  local = TRUE)$value
  source(file.path("server", "cell_annotation.server.R"),  local = TRUE)$value
  source(file.path("server", "modalities.server.R"),  local = TRUE)$value
  source(file.path("server", "custom_metadata.server.R"),  local = TRUE)$value
  source(file.path("server", "file_conversion.server.R"),  local = TRUE)$value
  source(file.path("server", "compare_annotations.server.R"),  local = TRUE)$value
  source(file.path("server", "scibet.server.R"),  local = TRUE)$value
  source(file.path("server", "differential_expression.server.R"),  local = TRUE)$value
  
} # server end
