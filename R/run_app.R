# load required packages 
pkgs <- c('shiny','Seurat','ggplot2','ggthemes','BiocManager','plotly','cowplot', 'dplyr','devtools','shinyjqui','shinythemes','gtools','pryr','shinyFiles')
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

# run app 
runApp(appDir = '/Users/jsjoyal/Desktop/SCAP/R/shiny/')
