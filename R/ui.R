#' Shiny app ui object
#'
#' @import shiny
#' @import shinycssloaders
#' @import shinyFiles
#' @import shinyjqui
#' @import shinythemes
ui <- navbarPage(
  fluid = TRUE,
  theme = shinytheme('cosmo'),
  title = "Single Cell Analysis Portal",
  tabPanel(
    "Main",
    value = 'main',
    wellPanel(
      fluidRow(
        column(
          width = 4,
          tags$head(
            tags$style(
              HTML(
                "#inputs-table {border-collapse: collapse;}
                #inputs-table td {padding: 10px;vertical-align: bottom;}"
              ) #/ HTML
            ) #/ style
          ), #/ head
          tags$table(
            id = 'inputs-table',
            style = "width: 100%",
            tags$tr(
              tags$td(
                style = 'width: 50%',
                h3(textOutput('data_used'))
              ),
              tags$td(
                style = 'width: 50%; text-align: left;',
                div(
                  class = 'from-group shiny-input-container',
                  shinyDirButton("directory", "Select Dataset", "Please select a Dataset")
                )
              )
            )
          )
        ),
        column(
          width = 4, 
          uiOutput('assay_1') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '100px')
        ),
        column(
          width = 4,
          uiOutput('grouping_1') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '100px')
        )
      )
    ),
    conditionalPanel(
      condition = 'input.assay_1',
      fluidRow(
        style='padding:50px;',
        column(
          width = 6,
          wellPanel(
            style  = 'background: white;',
            plotlyOutput('dimplot_1', height = '100%', width = '100%')%>% withSpinner(color="black"),
            uiOutput('reduction_1', style = 'padding: 10px')
          )
        ),
        column(
          width = 6,
          wellPanel(
            style  = 'background: white;',
            plotlyOutput('featureplot_1', height = '100%', width = '100%')%>% withSpinner(color="black"),
            uiOutput('featureplot_1_feature.select', style = 'padding: 10px')
          )
        ),
        div(style = "height:10px;")
      ),
      fluidRow(
        column(
          align = 'center',
          width = 12,
          offset = 0,
          wellPanel(
            style  = 'background: white;',
            fluidRow(
             uiOutput('dotplot_1') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '100px')
            ),
            fluidRow(
              column(
                width = 4,
                uiOutput('dotplot_1_feature_select', style = 'padding: 10px')
              ),
              column(
                width = 4,
                uiOutput('do_split')
              ),
              column(
                width = 4,
                conditionalPanel(
                  condition = "input.do_split=='yes'",
                  uiOutput('split_by'),
                  tags$head(tags$style("#split_by{overflow-y:scroll; max-height: 200px; background: white;}"))
                )
              )
            )
          )
        )
      )
    )
  ),
  tabPanel(
    "Cell Annotation",
    conditionalPanel(condition = '!input.assay_1', h2('Please Select Your Dataset on the Main Tab', style = 'text-align: center; font-style: italic;')),
    conditionalPanel(
      condition = 'input.assay_1',
      sidebarPanel(
        uiOutput('assay_2'),
        tabsetPanel(
          id = 'annot_panel',
          tabPanel(
            value = 'cell_annotation_cluster',
            "Annotate Clusters",
            uiOutput('grouping_2'),
            h2(textOutput('annotation_title')),
            uiOutput('annotations')
          ),
          tabPanel(
            value = 'cell_annotation_custom',
            "Custom Annotations",
            h2('Select groups of cells from either scatter plot and annotate them'),
            h3('Name selected cells'),
            fluidRow(
              column(
                width = 9,
                textInput('custom_name',label = NA, placeholder = "Enter cluster name here...")
              ),
              column(
                width = 3,
                actionButton("add_to_tmp", "Add Annotation", style='padding:4px; font-size:80%')
              )
            )
          )
        ),
        textInput('new_scheme_name', 'Name your annotation scheme', placeholder = "Enter scheme name here..."),
        actionButton("set_cell_types", "Rename Clusters/Cells"),
        h4("Selected Cells"),
        tableOutput('selected_cells'),
        tags$head(tags$style("#selected_cells{overflow-y:scroll; max-height: 200px; background: white;}"))
      ),
      mainPanel(
        fluidRow(
          column(
            width = 6,
            wellPanel(
              style  = 'background: white;padding: 20px',
              plotlyOutput('dimplot_2', height = '100%', width = '100%') %>% withSpinner(color="black"),
              uiOutput('reduction_2',style = 'padding:10px')
            )
          ),
          column(
            width = 6,
            wellPanel(
              style  = 'background: white;',
              plotlyOutput('featureplot_2', height = '100%', width = '100%') %>% withSpinner(color="black"),
              uiOutput('featureplot_2_feature_select',style = 'padding:10px')
            )
          )
        ),
        wellPanel(
          style  = 'background: white;',
          fluidRow(
            column(
              width = 8,
              fluidRow(
                h2('Find the defining markers of the selected cells', style = 'text-align: center;'),
                actionButton('find.markers', "Find Markers for Selected Cells"),
                h4('To clear selection... Double click on either plot and click button')
              ),
              fluidRow(
                dataTableOutput('markers') %>% withSpinner(color="black"),
                tags$head(tags$style("#markers{overflow-x:scroll; max-width: 90%; background: white;}"))
              )
            ),
            column(
              width = 4,
              h2("NCBI Gene Summary"),
              fluidRow(
                uiOutput('gene_query')
              ),
              fluidRow(
                actionButton('query_ncbi', label = 'Search')
                
              ),
              fluidRow(
                radioButtons('organism', label = 'Choose the organism', choices = c('Human', 'Mouse'))
              ),
              fluidRow(
                textOutput('gene_summary')
              )
            )
          )
        )
      )
    )
  ),
  tabPanel(
    "Custom Meta Data Groupings",
    conditionalPanel(condition = '!input.assay_1', h2('Please Select Your Dataset on the Main Tab', style = 'text-align: center; font-style: italic;')),
    conditionalPanel(
      condition = 'input.assay_1',
      sidebarPanel(
        tabsetPanel(
          id = "custom_meta_tab",
          tabPanel(
            value = "main",
            "Create Custom Groupings",
            uiOutput(outputId = 'assay_3'),
            uiOutput(outputId = 'meta_group_checkbox'),
            textInput(inputId = 'new_meta_group', placeholder = 'Enter meta data grouping name here...',label = NULL),
            actionButton(inputId = 'add_new_meta_group', label = "Add")
          ),
          tabPanel(
            value = "help",
            "Why would you want to do this?",
            h4("Sometimes it is useful to visualize the expression of genes within groups of cells other than common groupings, such as clusters, genotype, etc..."),
            h4("It may be especially useful to visualize gene expression within a combination of cell groups. For example, if we were originally provided meta data 
              that included the groupings for cell clusters (cluster_1, cluster_3, ...) and cell treatments (ctl, stim), we could visualize the gene expression
               between cell clusters or between treatments, but we may wish to also visualize the expression of a gene between treatments within each cluster."),
            h4("The potential usefulness of this can be visualized on the right. We see that Hepatocytes express Hrsp12. However, when we split up Hepatocytes 
               based on their treatment, we see that the stimulated Hepatocytes show greater expression of Hrsp12.")
          )
        )
      ),
      mainPanel(
        conditionalPanel(
          condition = "input.custom_meta_tab=='main'",
          h1(textOutput(outputId = 'example_meta_group_name')),
          tableOutput(outputId = 'example_meta_group')
        ),
        conditionalPanel(
          condition = "input.custom_meta_tab=='help'",
          fluidRow(
            column(
              width = 6,
              h4("Cell types"),
              img(src = "dotplot_celltype.png", height = "100%", width = "100%")
            ),
            column(
              width = 6,
              h4("Treatments"),
              img(src = "dotplot_treatment.png", height = "100%", width = "100%")
            )
          ),
          fluidRow(
            h4("Cell types + Treatments"),
            img(src = "dotplot_celltype_treatment.png", height = "50%", width = "50%")
          )
        )
      )
    )
  ),
  tabPanel(
    "File Conversion",
    sidebarPanel(
      tabsetPanel(
        tabPanel(
          'Seurat to Loom',
          fluidRow(
            style = 'padding:10px',
            shinyFilesButton(id = "object_file", label = "Choose Seurat Object to Convert", title = "Choose Seurat Object to Convert", multiple = F, style = "width:100%")
          ),
          fluidRow(
            style = 'padding:10px',
            shinyDirButton(id = "new_directory", label = "Select Directory", title = "Choose where to save converted loom files. The slected directory should be a new directory or empty", style = "width:100%")
          ),
          fluidRow(
            style = 'padding:10px',
            withBusyIndicatorUI(actionButton(inputId = 'sl_convert', label = 'Convert', icon = icon('arrow-circle-down'), style = "width:100%"))
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
          conditionalPanel(condition = '!input.assay_1', h2('Please Select Your Dataset on the Main Tab', style = 'text-align: center; font-style: italic;')),
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
        )
      )
    ),
    mainPanel(
      #h1('Stuff Goes Here')
    )
  )
)# end ui
#jqui_resizable(jqui_draggable(
#%>% withSpinner(color="#0dc5c1")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++#