#   ####
#   #-- Functions for annotation comparison--#
#   ####
#-- get the first list of potential annotations from loom--#
output$comp_anno_list1 <- renderUI({
  req(input$assay_1)
  req(rvalues$obs)
  req(rvalues$obs_cat)
  
  ## get annotation options from rvalues
  annotation_options <- rvalues$obs[rvalues$obs_cat]
  
  #group.by <- list(rvalues$h5ad[[1]]$obs[input$grouping_1][,,drop=TRUE])
  selectInput(
    inputId = 'comp_anno_1',
    label = 'Select the annotation you want to compare!',
    choices = annotation_options,
    multiple = FALSE)
})

#-- second list of annotations to compare against--#
output$comp_anno_list2 <- renderUI({
  req(input$assay_1,input$comp_anno_1)
  
  annotation_options <- rvalues$obs[rvalues$obs_cat]
  annotation_options <- annotation_options[!grepl(input$comp_anno_1,annotation_options)]
  
  selectInput(
    inputId = 'comp_anno_2',
    label = 'Select the annotation to compare against!',
    choices = annotation_options,
    multiple = FALSE)  
})


## Function to create data table for sankey diagram
sankey_comp <- eventReactive(input$compare_annos, {
  req(input$comp_anno_1, input$comp_anno_2)
  
  annos_to_compare <- rvalues$h5ad[[1]]$obs[c(input$comp_anno_1, input$comp_anno_2)]
  
  annos_to_compare_stats <- annos_to_compare %>%
    group_by(get(input$comp_anno_1),get(input$comp_anno_2)) %>%
    tally() %>%
    ungroup()
  
  colnames(annos_to_compare_stats) <- c("anno1","anno2","n")
  
  annos_to_compare_stats <- as.data.frame(annos_to_compare_stats)
  rownames(annos_to_compare_stats) <- 1:nrow(annos_to_compare_stats)
  
  # annos_to_compare_stats$anno1 <- as.character(annos_to_compare_stats$anno1)
  # annos_to_compare_stats$anno2 <- as.character(annos_to_compare_stats$anno2)
  
  annos_to_compare_stats <- annos_to_compare_stats %>%
    mutate("anno1" = paste(input$comp_anno_1,anno1,sep=": ")) %>%
    mutate("anno2" = paste(input$comp_anno_2,anno2,sep=": "))
  
  joined_annos_1 <- unique(annos_to_compare_stats$anno1) 
  joined_annos_2 <- unique(annos_to_compare_stats$anno2)
  
  annos_to_compare_stats$IDsource= match(annos_to_compare_stats$anno1, joined_annos_1) - 1
  annos_to_compare_stats$IDtarget= (match(annos_to_compare_stats$anno2, joined_annos_2))  - 1
  annos_to_compare_stats <- annos_to_compare_stats %>%
    mutate("IDtarget" = IDtarget + length(joined_annos_1))
  
  return(annos_to_compare_stats)
})


#### Compare annotations elements
## Sankey diagram to compare two annotations
output$sankey_diagram <- renderPlotly({
  req(sankey_comp())
  
  joined_annos <- c(sankey_comp()$anno1,sankey_comp()$anno2)
  joined_annos <- unique(joined_annos)
  
  p <- plot_ly(
    type = "sankey",
    orientation = "h",
    
    node = list(
      label = joined_annos,
      pad = 15,
      thickness = 40,
      line = list(
        color = "black",
        width = 0.5
      )
    ),
    
    link = list(
      source = sankey_comp()$IDsource,
      target = sankey_comp()$IDtarget,
      value =  sankey_comp()$n
    )
  ) %>%
    layout(font = list(
      size = 18)
    )
  
  p
  
}) # sankey_diagram end