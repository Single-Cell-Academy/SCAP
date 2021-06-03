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
        tags$head(
          tags$style(
            type="text/css",
            "#image img {max-width: 100%; width: auto; height: auto}"
          )
        ),
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
                shinyFilesButton("h5ad_in", "Select an H5ad Dataset", "Please select an H5ad Dataset", multiple = FALSE)
              )
            )
          )
        )
      ),
      ## Only show image when no data is loaded!
      conditionalPanel(
        condition = "!input.assay_1",	 
        column(
          width = 12, 
          align="center",
          imageOutput("genap_logo", height = "80%")
        )
      ),
      ## Dropdown menus on main page
      column(
        width = 4,
        uiOutput('assay_1') #%>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '100px')
      ),
      column(
        width = 4,
        uiOutput('grouping_1') #%>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '100px')
      )
    ) # End of this fluid row
  ),
  ## Panel for Embedding plot and Featureplot
  conditionalPanel(
    condition = 'input.assay_1',
    fluidRow(
      style='padding:50px;',
      column(
        width = 6,
        wellPanel(
          style  = 'background:white;',
          fluidRow(
            plotlyOutput('dimplot_1', height = '400px') %>% withSpinner(color="black", proxy.height = '400px')
          ),
          fluidRow(
            uiOutput('reduction_1', style = 'padding: 10px')
          )
        )
      ),
      column(
        width = 6,
        wellPanel(
          style  = 'background:white;',
          fluidRow(
            plotlyOutput('featureplot_1', height = '400px') %>% withSpinner(color="black", proxy.height = '400px')
          ),
          fluidRow(
            column(
              width = 3,
              br(),
              radioButtons('nebulosa_on', label = 'Nebulosa plot', choices = c('yes', 'no'), selected = 'no', inline = TRUE)
            ),
            column(
              width = 9,
              uiOutput('featureplot_1_feature_select', style = 'padding: 10px')
            )
          )
        )
      ),
      div(style = "height:10px;")
    ),
    conditionalPanel( ## Only show Dotplot when user selects a categorical variable
      condition = "output.grouping_1_type=='yes'",
      fluidRow(
        column(
          align = 'center',
          width = 12,
          offset = 0,
          wellPanel(
            style  = 'background: white;',
            fluidRow(
              uiOutput('dotplot_1') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '400px')
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
      ) # end of conditional panel for categorical variable
    )
  )
)