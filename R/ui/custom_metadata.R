## Custom metadata panel start
tabPanel(
  "Custom Metadata",
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
)