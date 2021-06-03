## Custom metadata panel end  
#### Scibet panel
tabPanel(
  "SciBet",
  conditionalPanel(
    condition = '!input.assay_1', 
    h2('Please Select Your Dataset on the Main Tab', 
       style = 'text-align: center; font-style: italic;')
  ),
  conditionalPanel(
    condition = 'input.assay_1',
    sidebarPanel(
      h3("Choose your reference dataset:",align="left"),
      selectInput(inputId="sci_bet_species", label="Choose a species!",
                  choices=c("Mouse","Human"),
                  selected="Mouse", multiple=FALSE, selectize=FALSE),
      uiOutput("scibet_reference_sets"),
      uiOutput("predict_cells_button"),
      br(),
      h3("Save predictions:"),
      br(),
      uiOutput("add_predictions_button"),
    ),
    # Show a plot of the generated distribution
    mainPanel(
      # # Info box containing information about the selected dataset
      h1("Explanation"),
      p("You can use this Panel to predict cell types in your data from published, annotated Datasets using SciBet. See the table below,
        for more information on available references:"),
      tags$a(href="https://www.nature.com/articles/s41467-020-15523-2",target = "blank", "SciBet Manuscript!"),
      h2("References available"),
      reactableOutput('scibet_references', height = 500),
      h1("Results"),
      h2("Distribution of predicted cell types:"),
      plotOutput("scibet_predictions_plot"),
    )
  )
)