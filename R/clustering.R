library(tidyverse)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(plotly)
library(shinybusy)

clusteringTabServer <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    # Reactive values to store data
    selected_variable_genes <- reactiveVal()
    samples_anova_results <- reactiveVal()
    
    # Apply filters based on user inputs
    observeEvent(input$run_analysis, {
      isolate({
        # Only recalculates when the above reactive values change
        samples_anova_results(calculate_variable_genes())
      })

    })

    calculate_variable_genes <- function(){
      req(dataset$filtered_data(), dataset$groups_data)

      df <- dataset$filtered_data()
      groups <- dataset$groups_data()
      
      # remove this column for simplicity 
      gene_data <- df %>% select(-Symbols)
      colnames(gene_data)[1] <- "Genes"
      colnames(groups)[1] <- "Condition"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"
      
      
      gene_data_long <- gene_data %>%
        pivot_longer(
          cols = -Genes,  # Exclude Gene_Symbol from the melting process
          names_to = "Condition",  # This will contain the sample/condition names
          values_to = "Expression"  # This will contain the expression values
        )
      
      # Join the dataframes to map Exp_Grp values to Conditions in gene_data_long
      gene_data_long <- gene_data_long %>%
        left_join(groups, by = "Condition")
      
      group_means <- gene_data_long %>%
        group_by(Genes, EXPERIMENTAL_GROUP) %>%
        summarise(MeanExpression = mean(Expression, na.rm = TRUE), .groups = 'drop')
      
      max_group_means <- group_means %>%
        group_by(Genes) %>%
        summarise(MaxMean = max(MeanExpression))
      
      
      anova_results <- gene_data_long %>%
        group_by(Genes) %>%
        do({
          model <- lm(Expression ~ EXPERIMENTAL_GROUP, data = .)
          anova_model <- anova(model)
          data.frame(p_value = anova_model$`Pr(>F)`[1])
        }) %>%
        ungroup()
      
      anova_results$p_value_adj <- p.adjust(anova_results$p_value, method = "BH")
      anova_results %>%
        left_join(max_group_means, by = "Genes")
      
    }
    
    output$variable_genes_plot <- renderPlot({
      req(samples_anova_results())
      result_data <- samples_anova_results()
    
      # update selected variable genes
      variable_genes <- filter(result_data, p_value_adj < input$y_cutoff, MaxMean > input$x_cutoff)
      selected_variable_genes(variable_genes)
      
      # highlight selected variable genes in the plot
      result_data$highlight <- ifelse(result_data$p_value_adj < input$y_cutoff & result_data$MaxMean > input$x_cutoff, "Selected", "Unselected")
      
      result_data$MaxMean_log2 <- log2(result_data$MaxMean + 1)
      
      p <- result_data %>%
        ggplot(aes(x = MaxMean_log2, y = p_value_adj, label = Genes, color = highlight, text = paste("Genes:", Genes, "<br>x:", MaxMean, "<br>y:", p_value))) + 
        geom_point(size=3) +
        scale_color_manual(values = c("Selected" = "blue", "Unselected" = "grey")) +
        labs(x=paste0("Max Mean"),
             y=paste0("P Values Adjusted")) +
        theme_minimal() + 
        theme(legend.position = "right") +
        theme_linedraw(base_size=16) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              legend.title = element_blank(),
              legend.text = element_text(size = 10),
              legend.position = "right")
      p
      # ggplotly(p, tooltip = "text")
    })
    
    output$number_of_genes <- renderText(paste0("Number of Selescted Genes: ", nrow(selected_variable_genes())))
    
    create_cluster_plot <- function(data, cluster_options, max_clusters) {
      data_sample_expr <- data[,-c(1)]
      gene_expression_z <- scale(data_sample_expr)
      gene_expression_z[is.na(gene_expression_z)] <- 0 
      rownames(gene_expression_z) <- data$Genes

      if (cluster_options == 1) {
        set.seed(40)
        clusters <- kmeans(gene_expression_z, centers=max_clusters, nstart=25)
        
        annotation_df <- data.frame(Genes = data$Genes, Clusters = clusters[["cluster"]])
        rownames(annotation_df) <- data$Genes
        annotation_df <- annotation_df %>% arrange(Clusters) %>% select(-Genes)
        
        order <- order(clusters$cluster)
        cor_matrix_ordered <- gene_expression_z[order, ]
        p <- pheatmap::pheatmap(cor_matrix_ordered, annotation_row = annotation_df,
                                                               cluster_rows = FALSE,
                                                               cluster_cols = FALSE,
                                                               fontsize_row = 10,
                                                               main = "K-means Clustered Matrix",
                                                               treeheight_row = 30,
                                                               # labels_row = FALSE,
                                                               angle_col = 315,
                                                               fontsize_col = 10,
                                                               width = 10,
                                                               height = 10)
        
        return(p)
      }
      
      p <- pheatmap::pheatmap(gene_expression_z, 
                             cluster_rows = TRUE,
                             cluster_cols = TRUE, 
                             clustering_distance_rows = "correlation",
                             clustering_distance_cols = "correlation",
                             fontsize_row = 10,
                             main = "Hierarchically Clustered Matrix",
                             treeheight_row = 30,
                             angle_col = 315,
                             fontsize_col = 10,
                             width = 10,
                             height = 10)
      return(p)
    }
    
    output$heatmap_plot <- renderPlot({
      req(dataset$filtered_data(), selected_variable_genes())
      variable_genes <- selected_variable_genes()
      df <- dataset$filtered_data()
      gene_data <- df %>% select(-Symbols)
      heatmap_df <- gene_data %>%
        inner_join(variable_genes, by = "Genes") %>%
        select(-p_value, -p_value_adj, -MaxMean)
      
      p <- create_cluster_plot(heatmap_df, input$cluster_options, input$clustering_k)
      print(p)
      return(p)
    })
    
    
  })
}

clusteringTabUI <- function(id) {
  fluidPage(
    titlePanel("Clustering of Gene Expression"),
    fluidRow(
      column(2, actionButton(NS(id, "run_analysis"), "Run Analysis"))
    ),
    fluidRow(
      column(2, numericInput(NS(id, "x_cutoff"), "Max Mean Cutoff", value = 0)),
      column(2, numericInput(NS(id, "y_cutoff"), "P Values Adjusted Cutoff", value = 1)),
      column(2, selectInput(NS(id, "cluster_options"), "Clustering Methods",
                            choices = c("k means" = 1, "hierarchical" = 2), selected = 1)),
      column(2, numericInput(NS(id, "clustering_k"), "K Values", value = 5))
    ),
    fluidRow(
      column(6,
             textOutput(NS(id, "number_of_genes")),
             add_busy_spinner(spin = "fading-circle", color = "#000000"),
             plotOutput(NS(id, "variable_genes_plot"), width = "100%", height = "700px")),
      column(6, plotOutput(NS(id, "heatmap_plot"), width = "100%", height = "700px"))
    ),
  )
  
}
