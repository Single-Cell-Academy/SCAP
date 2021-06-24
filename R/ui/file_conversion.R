### File conversion panel
tabPanel(
  "File Conversion",
  sidebarPanel(
    fluidRow(
      style = 'padding:10px',
      column(
        width = 6,
        selectInput('from', label = 'From', choices = c('rds', 'h5ad', 'loom'), selected = 'rds', multiple = FALSE)
      ),
      column(
        width = 6, 
        selectInput('to', label = 'To', choices = c('rds', 'h5ad'), 'h5ad', multiple = FALSE)
      )
    ),
    fluidRow(
      style = 'padding:10px',
      shinyFilesButton(id = "file_1", label = "Choose a Seurat (.rds), Scanpy (.h5ad) or Loom (.loom) Object to Convert", title = "Choose a Seurat (.rds), Scanpy (.h5ad) or Loom (.loom) Object to Convert", multiple = F, style = "width:100%")
    ),
    fluidRow(
      style = 'padding:10px',
      shinyDirButton(id = "out_dir_2", label = "Select Directory", title = "Choose where to save the converted file", style = "width:100%")
    ),
    fluidRow(
      style = 'padding:10px',
      textInput('out_file_2', label = 'Name converted file.\nIf multiple assays are present, assay names will automatically be appended to chosen file name')
    ),
    fluidRow(
      style = 'padding:10px',
      actionButton(inputId = 'scap_convert', label = 'Convert', icon = icon('arrow-circle-down'), style = "width:100%")
    )
  ),
  mainPanel(
    fluidRow(
      conditionalPanel(
        condition = 'input.file_1',
        h3('Selected input file:'),
        h4(textOutput(outputId = 'chosen_in_file'))
      )
    ),
    fluidRow(
      conditionalPanel(
        condition = 'input.out_dir_2',
        h3('Selected output directory:'),
        h4(textOutput(outputId = 'chosen_out_directory'))
      )
    ),
    fluidRow(
      conditionalPanel(
        condition = 'input.out_file_2',
        h3('Selected output file:'),
        h4(textOutput(outputId = 'chosen_out_file'))
      )
    )
  )
)