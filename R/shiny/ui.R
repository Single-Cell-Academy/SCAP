# User Interface

ui <- navbarPage(
  fluid = TRUE,
  theme = shinytheme('cosmo'),
  title = "Single Cell Analysis Portal",
  tabPanel(
    "Main",
    sidebarPanel(
      fluidRow(
        h3(textOutput('data_used')),
        shinyDirButton("directory", "Select Dataset", "Please select a Dataset")
      ),
      fluidRow(
        uiOutput('assay_1') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '100px'),
        uiOutput('reduction_1') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '100px'),
        uiOutput('grouping_1') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '100px'),
        uiOutput('featureplot_1_feature.select') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '100px'),
        uiOutput('dotplot_1_feature_select') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '100px')
      )
    ),
    mainPanel(
      "",
      fluidRow(style='padding:100px;',
               column(
                 width = 6,
                 plotlyOutput('dimplot_1', height = '100%', width = '100%')%>% withSpinner(color="black")
               ),
               column(
                 width = 6,
                 plotlyOutput('featureplot_1', height = '100%', width = '100%')%>% withSpinner(color="black")
               ),
               div(style = "height:10px;")
      ),
      fluidRow(
        column(
          width = 12,
          offset = 0,
          plotOutput('dotplot_1', height = '500px', width = '100%')%>% withSpinner(color="black")
        )
      )
    )
  ),
  tabPanel(
    "Cell Annotation",
    sidebarPanel(
      uiOutput('assay_2'),
      tabsetPanel(
        id = 'annot_panel',
        tabPanel(
          value = 'cluster',
          "Annotate Clusters",
          uiOutput('grouping_2'),
          #uiOutput('annots'),
          textInput('new_name', 'Name your annotation scheme', placeholder = "Enter scheme name here..."),
          actionButton("setCellsTypes", "Rename Clusters/Cells")
        ),
        tabPanel(
          value = 'custom',
          "Custom Annotations",
          textInput('custom_name','Name Selected Cluster', placeholder = "Enter cluster name here..."),
          actionButton("add_to_tmp", "Add Annotation"),
          textInput('custom_scheme','Name Custom Scheme', placeholder = "Enter scheme name here..."),
          actionButton("save_custom_scheme", "Save Custom Annotations")
        )
      ),
      uiOutput('reduction_2'),
      uiOutput('featureplot_2_feature_select'),
      h4("Selected Cells"),
      tableOutput('selected_cells'),
      tags$head(tags$style("#cells{overflow-y:scroll; max-height: 200px; background: ghostwhite;}"))
    ),
    mainPanel(
      fluidRow(
        style='padding:100px;',
        column(
          width = 6,
          plotlyOutput('dimplot_2', height = '100%', width = '100%')%>% withSpinner(color="black")
        ),
        column(
          width = 6,
          plotlyOutput('featureplot_2', height = '100%', width = '100%')%>% withSpinner(color="black")
        )
      ),
      fluidRow(
        style='padding:50px;',
        actionButton('find.markers', "Find Markers for Selected Cells"),
        h4('To clear selection... Double click on either plot and click button')
      ),
      fluidRow(
        tableOutput('markers') %>% withSpinner(color="black")
      )
    )
  )
  #,
  # tabPanel(
  #   "Create Custom Meta Data Groupings",
  #   sidebarPanel(
  #     uiOutput(outputId = 'assay.3'),
  #     uiOutput(outputId = 'meta_group_checkbox'),
  #     textInput(inputId = 'new_meta_group', placeholder = 'Enter meta data grouping name here...',label = NULL),
  #     actionButton(inputId = 'add_new_meta_group', label = "Add")
  #   ),
  #   mainPanel(
  #     h1(textOutput(outputId = 'example_meta_group_name')),
  #     textOutput(outputId = 'example_meta_group')
  #   )
  # )
)
#jqui_resizable(jqui_draggable(
#%>% withSpinner(color="#0dc5c1")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++#