## Custom metadata panel start
tabPanel(
  "Differential expression",
  conditionalPanel(condition = '!input.assay_1', h2('Please Select Your Dataset on the Main Tab', style = 'text-align: center; font-style: italic;')),
  conditionalPanel(
    condition = 'input.assay_1',
    sidebarPanel(
      tabsetPanel(
        tabPanel(
          value = "main",
          "DE testing",
          br(),
          uiOutput(outputId = 'de_annotation_list'),
          uiOutput(outputId = 'de_group_1_list'),
          uiOutput(outputId = 'de_group_2_list'),
          actionButton(inputId = 'run_de_analysis', label = "Run differential expression!")),
        tabPanel(
          value = "help",
          "What can I do on this page?",
          br(),
          p("On this page you can perform a simple differential expression analysis between two groups of cells. First, chose the observations
               that you want to compare. If you require grouping of multiple observations for comparison, you can add these via custom metadata. Then,
               select the two groups of cells you would like to compare and click run differential expression."),
          br(),
          p("The differential expression analysis implemented here will use the wilcoxauc() function from the presto package. Please note, that while
            this approach generally reveals strong effects, this is a very simple way to perform differential expression analysis, that does not take into account any information on technical or 
            biological replicates and can therefore have inflated statistics. To run more precise differential expression models, please consider
              approaches outside of SCAP.")
        )
      )
    ),
    mainPanel(
      ## Message shown at the very beginning

      shinyjs::hidden(
        div(id = "empty_de",
          wellPanel(h2('Please select your settings for differential expression analysis on the left.', 
             style = 'text-align: center; font-style: italic;'))
        )),
      
      shinyjs::hidden(
        div(
          id = "de_results",
                       fluidRow(
                         column(width = 8,
                                reactableOutput("de_res_table"),
                         ),
                         column(width = 4,
                                plotOutput("de_violin_plot"))
                       ),
                       br(),
                       br(),
                       fluidRow(
                         column(width = 6,
                                plotOutput("de_volcano_plot")
                         ),
                         column(width = 6,
                                plotOutput("de_avg_exp_plot"))
                       )
          )
        )
    )
  )
)