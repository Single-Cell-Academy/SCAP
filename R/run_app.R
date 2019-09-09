#' Run SCAP
#' @export
SCAP <- function() {
  shinyApp(ui = ui, server = server)
}



#-- Install and load required packages --#
# pkgs <- c('shiny','shinyjs','Seurat','ggplot2','ggthemes','BiocManager','plotly','cowplot', 'dplyr','devtools','shinyjqui','shinythemes','gtools','pryr','shinyFiles','Matrix','hdf5r','MODIS')
# lapply(X = pkgs, FUN = function(x){
#   if(!require(x, character.only = TRUE))
#     install.packages(x, repos = 'https://cloud.r-project.org/')
# })
# lapply(X = pkgs, FUN = function(x){
#   if(!require(x, character.only = TRUE))
#     library(x)
# })
# if(!require('shinycssloaders'))
#   devtools::install_github('https://github.com/andrewsali/shinycssloaders')
# library(shinycssloaders)
# if(!require('presto'))
#   devtools::install_github('https://github.com/immunogenomics/presto.git')
# library(presto)
# library(shinycssloaders)
# if(!require('loomR'))
#   devtools::install_github('https://github.com/mojaveazure/loomR.git', ref = 'develop')
# library(loomR)
# 
# #-- Source required functions --#
# source("/Users/jsjoyal/Desktop/SCAP/R/SCAP_functions.R")

#runApp(appDir = '/Users/jsjoyal/Desktop/SCAP/R/shiny/')