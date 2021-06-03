tabPanel(
  "Modalities",
  value = 'modalities',
  conditionalPanel(
    condition = '!input.assay_1', 
    h2('Please Select Your Dataset on the Main Tab', style = 'text-align: center; font-style: italic;')
  ),
  conditionalPanel(
    condition = 'input.assay_1',
    wellPanel(
      fluidRow(
        column(
          width = 2,
          uiOutput('data_used_mod')
        ),
        column(
          width = 2,
          shinyFilesButton("h5ad_in_mod", "Select an H5ad Dataset", "Please select an H5ad Dataset", multiple = FALSE)
        ),
        column(
          width = 2,
          uiOutput('assay_mod') 
        ),
        column(
          width = 2,
          uiOutput('grouping_mod') 
        ),
        column(
          width = 2,
          selectInput(inputId = "mod_sel", label = "Select Modality Type", choices = c("None", "CRISPR", "CITE-Seq", "Feature Correlation"), selected = "None")
        )
      )
    ),
    conditionalPanel(
      condition = "input.assay_1",
      fluidRow(
        style='padding:50px;',
        column(
          width = 6,
          wellPanel(
            style  = 'background:white;',
            fluidRow(
              plotlyOutput('dimplot_mod', height = '400px') %>% withSpinner(color="black", proxy.height = '400px')
            ),
            fluidRow(
              uiOutput('reduction_mod', style = 'padding: 10px')
            )
          )
        ),
        column(
          width = 6,
          wellPanel(
            style  = 'background:white;',
            fluidRow(
              plotlyOutput('featureplot_mod', height = '400px') %>% withSpinner(color="black", proxy.height = '400px')
            ),
            fluidRow(
              column(
                width = 3,
                br(),
                radioButtons('nebulosa_mod_on', label = 'Nebulosa plot', choices = c('yes', 'no'), selected = 'no', inline = TRUE)
              ),
              column(
                width = 9,
                uiOutput('featureplot_mod_feature_select', style = 'padding: 10px')
              )
            )
          )
        ),
        div(style = "height:10px;")
      ),
      conditionalPanel( ## Show Ridgeplot when user selects CITE-seq
        condition = "input.mod_sel=='CITE-Seq'",
        fluidRow(
          column(
            align = 'center',
            width = 12,
            offset = 0,
            wellPanel(
              style  = 'background: white;',
              fluidRow(
                plotOutput('ridgeplot_mod') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '400px'),
                uiOutput('ridgeplot_mod_feature_select', style = 'padding: 10px')
              )
            )
          )
        )
      ),
      conditionalPanel( ## Show differential expression toolkit for CRISPR
        condition = "input.mod_sel=='CRISPR'",
        wellPanel(
          style  = 'background: white;',
          fixedRow(align = "center",
                   column(align = 'left',width = 2,
                          uiOutput('crispr_feature_1', style = 'padding: 10px'),
                   ),            
                   column(align = 'left',width = 2,
                          uiOutput('crispr_feature_1_slider', style = 'padding: 10px'),
                   ),          
                   column(align = 'left', width = 4,
                          plotOutput('crispr_feature_1_dist',
                                     height = '150px')
                   ),
                   column(align = 'left', width = 2,
                          textOutput("crispr_feature_1_cells_print")
                   )
          ),
          hr(),
          fixedRow(align = "center",
                   column(align = 'left',width = 2,
                          uiOutput('crispr_feature_2', style = 'padding: 10px'),
                   ),
                   column(align = 'left',width = 2,
                          uiOutput('crispr_feature_2_slider', style = 'padding: 10px'),
                   ), 
                   column(align = 'left', width = 4,
                          plotOutput('crispr_feature_2_dist',
                                     height = '150px')
                   ),
                   column(align = 'left', width = 2,
                          textOutput("crispr_feature_2_cells_print"),
                   ),
                   column(align = 'left', width = 2,
                          actionButton(inputId = 'crispr_de_analysis', label = "Compare RNA content of feature cells!"))
          )
        ),
        br(),
        shinyjs::hidden(
          wellPanel(id = "crispr_res",
                    style  = 'background: white;',
                    fluidRow(
                      column(align = "center", width = 4,
                             plotlyOutput("crispr_avg_gene_exp",
                                          height = 'auto') %>% withSpinner(type = 3,color.background = "white",
                                                                           color="black", proxy.height = '400px')),
                      
                      column(align = "center", width = 4,
                             reactableOutput("crispr_avg_gene_exp_tbl",
                                             height = 'auto') %>% withSpinner(type = 3,color.background = "white",
                                                                              color="black", proxy.height = '400px')),
                      
                      column(align = "center", width = 4,
                             plotOutput("crispr_gene_vlnplot") %>% withSpinner(type = 3,color.background = "white",
                                                                               color="black", proxy.height = '400px')
                      )
                    )
          )
        )
      ),
      conditionalPanel(
        condition = "input.mod_sel=='Feature Correlation'",
        sidebarPanel(
          uiOutput('corr_grouping'),
          uiOutput('corr_sub_grouping'),
          uiOutput('corr_fs_1'),
          uiOutput('corr_fs_2'),
          htmlOutput('corr_stats')
        ),
        mainPanel(
          plotOutput('corr_plot') %>% withSpinner(color="black", proxy.height = '400px')
        )
      )
    )
  )
)