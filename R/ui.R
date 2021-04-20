#' Shiny app ui object 

library("shiny")
library("shinycssloaders")
library("plotly")
library("reactable")
library("shinythemes")
library("shinyFiles")

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
  ),
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
  ),
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
            width = 4,
            tags$head(
              tags$style(
                HTML(
                  "#inputs-table-2 {border-collapse: collapse;}
                  #inputs-table-2 td {padding: 10px;vertical-align: bottom;}"
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
              id = 'inputs-table-2',
              style = "width: 100%",
              tags$tr(
                tags$td(
                  style = 'width: 50%',
                  h3(textOutput('data_used_mod'))
                ),
                tags$td(
                  style = 'width: 50%; text-align: left;',
                  div(
                    class = 'from-group shiny-input-container',
                    shinyFilesButton("h5ad_in_mod", "Select an H5ad Dataset", "Please select an H5ad Dataset", multiple = FALSE)
                  )
                )
              )
            )
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
            selectInput(inputId = "mod_sel", label = "Select Modality Type", choices = c("None", "CRISPR", "CITE-Seq"), selected = "None")
          )
        )
      ),
      conditionalPanel(
        condition = 'input.assay_mod',
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
        conditionalPanel( ## Show Ridgeplot when user selects CITE-seq and show differential expression toolkit for CRISPR
          condition = "input.mod_sel=='CITE-Seq'",
          fluidRow(
            column(
              align = 'center',
              width = 6,
              offset = 0,
              wellPanel(
                style  = 'background: white;',
                fluidRow(
                  plotOutput('ridgeplot_mod') %>% withSpinner(color="black",type = 7,size = 0.5,proxy.height = '400px')
                ),
                fluidRow(
                  column(
                    width = 4,
                    uiOutput('ridgeplot_mod_feature_select', style = 'padding: 10px')
                  )
                )
              )
            )
          )
        ),
        conditionalPanel( ## Show Ridgeplot when user selects CITE-seq and show differential expression toolkit for CRISPR
          condition = "input.mod_sel=='CRISPR'",
          fixedRow(
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
          fixedRow(
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
                   actionButton(inputId = 'crispr_de_analysis', label = "Compare features RNA content!")
            )
          ),
          br(),
          fixedRow(
            column(align = "left", width = 4,
                   plotlyOutput("crispr_avg_gene_exp",
                                height = "400px")),
            column(align = "left", width = 4,
                   reactableOutput("crispr_avg_gene_exp_tbl")),
            column(align = "left", width = 4,
                   tableOutput("selected")
                   # plotOutput("crispr_gene_vlnplot")
                   ),
          ),
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
  ),
  #### Compare annotations
  tabPanel(
    "Compare annotations",
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
  ),
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
  ),
  tabPanel(
  "News",
  includeMarkdown("../news/news.md")
  ) ## end news tabPanel
)   # end ui