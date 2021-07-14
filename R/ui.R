#' Shiny app ui object 

library("shiny")
library("shinycssloaders")
library("plotly")
library("reactable")
library("shinythemes")
library("shinyFiles")
library("shinyjs")
library("data.table")

version_number <-  strsplit(grep('Version', readLines("../DESCRIPTION"), fixed = TRUE, value = TRUE),split = " ")[[1]][2]

  ui <- navbarPage(

  fluid = TRUE,
  collapsible = TRUE,
  theme = shinytheme('cosmo'),
  title = paste("Single Cell Analysis Portal v.",version_number,sep =""),
  ## Dashboard components (TabPanels) found in /ui
  source(file.path("ui", "main.R"),  local = TRUE)$value,
  source(file.path("ui", "cell_annotation.R"),  local = TRUE)$value,
  source(file.path("ui", "modalities.R"),  local = TRUE)$value,
  source(file.path("ui", "custom_metadata.R"),  local = TRUE)$value,
  source(file.path("ui", "differential_expression.R"),  local = TRUE)$value,
  source(file.path("ui", "scibet.R"),  local = TRUE)$value,
  source(file.path("ui", "compare_annotations.R"),  local = TRUE)$value,
  source(file.path("ui", "file_conversion.R"),  local = TRUE)$value,
  source(file.path("ui", "changelog.R"),  local = TRUE)$value,
  
  useShinyjs()
)   # end ui