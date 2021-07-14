### File conversion panel
tabPanel(
  "File Conversion",
  sidebarPanel(
    tabsetPanel(
      tabPanel(
        'Seurat/Scanpy to Loom',
        fluidRow(
          style = 'padding:10px',
          shinyFilesButton(id = "object_file", label = "Choose Seurat (.rds) / Scanpy (/.h5/.h5ad) Object to Convert ", title = "Choose Seurat (.rds) / Scanpy (/.h5/.h5ad) Object to Convert", multiple = F, style = "width:100%")
        ),
        fluidRow(
          style = 'padding:10px',
          shinyDirButton(id = "new_directory", label = "Select Directory", title = "Choose where to save converted loom files. The slected directory should be a new directory or empty", style = "width:100%")
        ),
        fluidRow(
          style = 'padding:10px',
          actionButton(inputId = 'sl_convert', label = 'Convert', icon = icon('arrow-circle-down'), style = "width:100%")
        ),
        fluidRow(
          conditionalPanel(
            condition = 'input.object_file',
            h3('Selected file:'),
            h4(textOutput(outputId = 'chosen_object_file')),
            h4(textOutput(outputId = 'notes'))
          )
        ),
        fluidRow(
          conditionalPanel(
            condition = 'input.new_directory',
            h3('Selected directory:'),
            h4(textOutput(outputId = 'chosen_new_directory'))
          )
        )
      ),
      tabPanel(
        "Loom to Seurat",
        conditionalPanel(
          condition = '!input.assay_1', 
          h2('Please Select Your Dataset on the Main Tab', style = 'text-align: center; font-style: italic;')
        ),
        conditionalPanel(
          condition = 'input.assay_1',
          fluidRow(
            style = "padding:10px",
            h3("Update the meta data of a Seurat object or create a new one"),
            radioButtons('overwrite', label = "Overwrite original Seurat object with new meta data", choices = c("yes", "no"), inline = T),
            conditionalPanel(
              condition = "input.overwrite=='no'",
              textInput('new_file_name', label = 'File name', value = NA, placeholder = "Enter the name of the updated/new Seurat file")
            )
          ),
          fluidRow(
            style = "padding:10px",
            shinyFilesButton(id = "update_file", label = "Choose Original Seurat Object", title = "Choose Seurat Object to Convert", multiple = F, style = "width:100%")
          ),
          fluidRow(
            style = "padding:10px",
            shinyDirButton(id = "new_directory_2", label = "Select Directory", title = "Choose where to save the updated Seurat file", style = "width:100%")
          ),
          fluidRow(
            style = "padding:10px",
            actionButton(inputId = 'ls_convert', label = 'Convert', icon = icon('arrow-circle-down'), style = "width:100%")
          ),
          fluidRow(
            h3('Selected file:'),
            h4(textOutput(outputId = 'chosen_update_file')),
            h3('Selected directory:'),
            h4(textOutput(outputId = 'chosen_new_directory_2'))
          )
        )
      ),
      tabPanel(
        "To h5ad",
        fluidRow(
          style = "padding:10px",
          h3("Convert loom or Seurat to h5ad"),
          shinyFilesButton(id = "in_file_to_h5ad", label = "Choose File to Convert", title = "Choose File to Convert", multiple = F, style = "width:100%"),
          textInput('out_file_to_h5ad', label = 'File name for h5ad Output', value = NA, placeholder = "Save h5ad File As"),
          shinyDirButton(id = "out_dir_to_h5ad", label = "Select where to save h5ad file", title = "Select where to save h5ad file", style = "width:100%"),
          radioButtons('is_legacy', label = "Legacy Conversion?", choices = c("yes", "no"), selected = "no", inline = T),
          conditionalPanel(
            condition = "input.is_legacy=='yes'",
            shinyFilesButton(id = "legacy_file_to_h5ad", label = "Choose Original Seurat Object", title = "Choose Original Seurat Object for Legacy Conversion", multiple = F, style = "width:100%")
          ),
          actionButton(inputId = 'to_h5ad', label = 'Convert', icon = icon('arrow-circle-down'), style = "width:100%")
        ),
        fluidRow(
          h3('Selected file:'),
          h4(textOutput(outputId = 'chosen_in_file_to_h5ad')),
          h3('Selected out file:'),
          h4(textOutput(outputId = 'chosen_out_file_to_h5ad')),
          h3('Selected directory:'),
          h4(textOutput(outputId = 'chosen_out_dir_to_h5ad')),
          h3('Selected legacy file:'),
          h4(textOutput(outputId = 'chosen_legacy_file_to_h5ad'))
        )
      )
    )
  ),
  mainPanel(
    #h1('Stuff Goes Here')
  )
)