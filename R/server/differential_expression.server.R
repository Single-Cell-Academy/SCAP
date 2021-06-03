## Differential expression server elements
output$de_annotation_list <- renderUI({
  req(input$assay_1)
  req(rvalues$obs)
  req(rvalues$obs_cat)
  
  ## get annotation options from rvalues
  annotation_options <- rvalues$obs[rvalues$obs_cat]
  
  #group.by <- list(rvalues$h5ad[[1]]$obs[input$grouping_1][,,drop=TRUE])
  selectInput(
    inputId = 'de_anno_sel',
    label = 'Select the observations you want to compare!',
    choices = annotation_options,
    multiple = FALSE)
})

## Baseline group for differential expression comparison
output$de_group_1_list <- renderUI({
  req(input$assay_1,input$de_anno_sel)
  
  anno_choices <- unique(rvalues$h5ad[[1]]$obs[c(input$de_anno_sel)][,1])
  
  selectInput(
    inputId = 'de_group_1_sel',
    label = 'Select the annotation to use as baseline!',
    choices = anno_choices,
    multiple = FALSE)  
})

## Group to compare against for differential expression comparison
output$de_group_2_list <- renderUI({
  req(input$assay_1,input$de_anno_sel, input$de_group_1_sel)
  
  anno_choices_2 <- unique(rvalues$h5ad[[1]]$obs[c(input$de_anno_sel)][,1])
  anno_choices_2 <- setdiff(anno_choices_2,input$de_group_1_sel)
  
  selectInput(
    inputId = 'de_group_2_sel',
    label = 'Select the annotation to compare against!',
    choices = anno_choices_2,
    multiple = FALSE)  
})

## Data frame containin expression for cells of DE group1
de_analysis_group1 <- reactive({
  df_annos <- rvalues$h5ad[[1]]$obs[c(input$de_anno_sel)][,1]
  df_annos <- data.frame("anno" = df_annos)
  df_annos$cell_ids <- 1:nrow(df_annos)
  df_cells <- subset(df_annos,anno == input$de_group_1_sel)
  if(rvalues$raw_dtype == "NULL" | rvalues$raw_dtype == "counts"){ ## Check if raw exists
    df_cells_exp_anno <- rvalues$h5ad[[1]]$X[df_cells$cell_ids,]
  }else if(rvalues$raw_dtype == "normalized"){ ## Only if the other two conditions fail, use raw values to calculate differential expression
    df_cells_exp_anno <- rvalues$h5ad[[1]]$raw$X[df_cells$cell_ids,] 
  }
  df_cells_exp_anno
})

## Data frame containing expression for cells of DE group2
de_analysis_group2 <- reactive({
  df_annos <- rvalues$h5ad[[1]]$obs[c(input$de_anno_sel)][,1]
  df_annos <- data.frame("anno" = df_annos)
  df_annos$cell_ids <- 1:nrow(df_annos)
  df_cells <- subset(df_annos,anno == input$de_group_2_sel)
  if(rvalues$raw_dtype == "NULL" | rvalues$raw_dtype == "counts"){ ## Check if raw exists
    df_cells_exp_anno <- rvalues$h5ad[[1]]$X[df_cells$cell_ids,]
  }else if(rvalues$raw_dtype == "normalized"){ ## Only if the other two conditions fail, use raw values to calculate differential expression
    df_cells_exp_anno <- rvalues$h5ad[[1]]$raw$X[df_cells$cell_ids,] 
  }
  df_cells_exp_anno
})

## Event when user clicks button to compare differential expression
observeEvent(input$run_de_analysis,{
  
  de_res <- reactive({
    de_group1 <- isolate(de_analysis_group1())
    de_group2 <- isolate(de_analysis_group2())
    
    ## Merge both expression tables
    matrix_tr <- t(rbind(de_group1,de_group2))
    rownames(matrix_tr) <- rvalues$features
    
    ## Create a feature vector for both groups
    feature_vec_1 <- replicate(nrow(de_group1),isolate(input$de_group_1_sel))
    feature_vec_2 <- replicate(nrow(de_group2),isolate(input$de_group_2_sel))
    feature_vec <- c(feature_vec_1,feature_vec_2)
    
    de_res <- wilcoxauc(matrix_tr, feature_vec)
    de_res <- de_res %>%
      arrange(desc(auc)) %>%
      subset(group == isolate(input$de_group_2_sel)) %>%
      dplyr::select(-c(pct_in,pct_out,group,statistic))
    
    de_res
  })
  
  avg_exp <- reactive({
    req(de_res())
    de_group1 <- isolate(de_analysis_group1())
    de_group2 <- isolate(de_analysis_group2())
    avg_exp_group1 <- colMeans(de_group1)
    avg_exp_group2 <- colMeans(de_group2)
    avg_exp_df <- data.frame("group1" = avg_exp_group1,
                             "group2" = avg_exp_group2,
                             "gene" = rvalues$features)
    
    avg_exp_df <- avg_exp_df %>%
      mutate("significant" = if_else(gene %in% subset(de_res(),padj < 0.05)$feature,"yes","no"))
    avg_exp_df
  })
  
  output$de_res_table <- renderReactable({
    de_res_tbl_view <- de_res() 
    
    reactable(isolate(de_res_tbl_view),
              sortable = TRUE,
              searchable = TRUE,
              selection = "single")
  })
  
  selected_de <- reactive(getReactableState("de_res_table", "selected"))
  selected_de_gene <- reactive({
    isolate(de_res())[selected_de(),]$feature
  })
  
  ## Violin plot for differential expression
  output$de_volcano_plot <- renderPlot({
    req(de_res())
    de_res_tbl <- de_res() %>%
      mutate("padj" = if_else(padj == 0,.Machine$double.xmin,padj))
    
    volcano_plot <- ggplot(isolate(de_res_tbl),aes(logFC,-log10(padj))) +
      geom_point() +
      bbc_style() +
      theme(axis.title = element_text(size = 16, face = "bold")) +
      labs(title = "Volcano plot",
           x = "log2FC",
           y = "-log10(padj)") +
      scale_y_continuous(limits = c(0,max(-log10(de_res_tbl$padj))+ 30))
    if(!is.null(selected_de())){
      volcano_plot <- volcano_plot + 
        geom_point(data = subset(isolate(de_res_tbl), feature == selected_de_gene()),
                   size = 5, fill = "red",color=  "black",pch = 21)
    }
    volcano_plot
  })
  
  ## Correlation plot for differential expression
  output$de_avg_exp_plot <- renderPlot({
    group1_name <- isolate(input$de_group_1_sel)
    group2_name <- isolate(input$de_group_2_sel)
    avg_exp_plot <- ggplot(isolate(avg_exp()),aes(group1,group2)) +
      geom_point() +
      geom_abline(linetype = 2) +
      labs(title = "Average expression",
           x = group1_name,
           y = group2_name) +
      bbc_style() +
      theme(axis.title = element_text(size = 16, face = "bold"))
    
    if(!is.null(selected_de())){
      avg_exp_plot <- avg_exp_plot + 
        geom_point(data = subset(isolate(avg_exp()), gene == selected_de_gene()),
                   size = 5, fill = "red",color=  "black",pch = 21)
    }
    avg_exp_plot
  })
  
  de_violin_data <- reactive({
    
    ## Merge both expression tables
    de_group_1 <- isolate(de_analysis_group1())
    group1_name <- isolate(input$de_group_1_sel)
    de_group_2 <- isolate(de_analysis_group2())
    group2_name <- isolate(input$de_group_2_sel)
    
    matrix <- rbind(de_group_1,de_group_2)
    colnames(matrix) <- rvalues$features
    matrix_sub <- matrix[,selected_de_gene()]
    
    ## Create a feature vector for both groups
    feature_vec_1 <- replicate(nrow(de_group_1),group1_name)
    feature_vec_2 <- replicate(nrow(de_group_2),group2_name)
    feature_vec <- c(feature_vec_1,feature_vec_2)
    
    final_df <- data.frame("exp" = matrix_sub,
                           "group" = feature_vec)
    
    ## set order of the two groups, such that group1 is always first
    final_df$group <- factor(final_df$group,
                             levels = c(group1_name,group2_name))
    
    final_df
  })
  
  
  ## Correlation plot for differential expression
  output$de_violin_plot <- renderPlot({
    req(selected_de_gene())
    ggplot(isolate(de_violin_data()),aes(group,exp, fill = group)) +
      geom_violin() +
      stat_summary(fun=mean, geom="point", size=5, color = "black") +
      scale_fill_manual(values = c("#4682B4","#B47846")) +
      bbc_style() +
      theme(axis.title = element_text(size = 16, face = "bold")) +
      theme(legend.position = "none") +
      labs(x = "Feature",
           y = "Gene expression",
           title = selected_de_gene())
  })
})