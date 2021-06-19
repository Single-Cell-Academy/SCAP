#### Compare annotations
tabPanel(
  "Compare Annotations",
  conditionalPanel(
    condition = '!input.assay_1', 
    h2('Please Select Your Dataset on the Main Tab', 
       style = 'text-align: center; font-style: italic;')
  ),
  conditionalPanel(
    condition = "input.assay_1",  
    fluidRow(
      column(
        width = 4,
        uiOutput("comp_anno_list1")
      ),
      column(
        width = 4,
        uiOutput("comp_anno_list2")
      ),
      column(
        width = 4,
        br(),
        actionButton("compare_annos","Compare annotations!")
      )
    ),
    fluidRow(
      column(
        width = 12,
        plotlyOutput("sankey_diagram", height = "auto")
      )
    )
  )
)