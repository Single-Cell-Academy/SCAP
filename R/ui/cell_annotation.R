tabPanel(
  "Cell Annotation",
  conditionalPanel(
    condition = '!input.assay_1', 
    h2('Please Select Your Dataset on the Main Tab', style = 'text-align: center; font-style: italic;')
  ),
  conditionalPanel(
    condition = 'input.assay_1',
    fluidRow(
      sidebarPanel(
        style = "overflow-y:scroll; max-height: 545px; min-height: 545px; background: white;",
        uiOutput('assay_2'), ###! replace by assay_1
        tabsetPanel(
          id = 'annot_panel',
          tabPanel(
            value = 'cell_annotation_cluster',
            "Annotate Clusters",
            uiOutput('grouping_2'),
            h2(textOutput('annotation_title')),
            uiOutput('annotations'),
            tags$head(tags$style("#annotations{overflow-y:scroll; max-height: 200px; background: white; border-color: black;}"))
          ),
          tabPanel(
            value = 'cell_annotation_custom',
            "Custom Annotations",
            h3('Select groups of cells from either scatter plot and annotate them'),
            h3('Name selected cells'),
            fluidRow(
              column(
                width = 9,
                textInput('custom_name',label = NULL, placeholder = "Enter cluster name here...")
              ),
              column(
                width = 3,
                actionButton("add_to_tmp", "Add Annotation", style='padding:4px; font-size:80%')
              )
            )
          )
        ),
        br(),
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
              plotlyOutput('dimplot_2', height = '400px') %>% withSpinner(color="black", proxy.height = '400px'),
              uiOutput('reduction_2',style = 'padding:10px')
            )
          ),
          column(
            width = 6,
            wellPanel(
              style  = 'background: white;',
              plotlyOutput('featureplot_2', height = '400px') %>% withSpinner(color="black", proxy.height = '400px'),
              uiOutput('featureplot_2_feature_select',style = 'padding:10px')
            )
          )
        )
      )
    ),
    fluidRow(
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
              DT::dataTableOutput('markers') %>% withSpinner(color="black"),
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
)