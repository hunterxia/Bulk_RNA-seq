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

      return(anova_results)
    }
    
    output$variable_genes_plot <- renderPlot({
      req(samples_anova_results())
      anova_results <- samples_anova_results()
      anova_results$p_value_adj <- p.adjust(anova_results$p_value, method = "BH")
      result_data <- anova_results %>%
        left_join(max_group_means, by = "Genes")
      
      # update selected variable genes
      variable_genes <- filter(result_data, p_value_adj > input$y_cutoff, MaxMean > input$x_cutoff)
      selected_variable_genes(variable_genes)
      
      # highlight selected variable genes in the plot
      result_data$highlight <- ifelse(result_data$p_value_adj > input$y_cutoff & result_data$MaxMean > input$x_cutoff, "Highlighted", "Normal")
      
      p <- result_data %>%
        ggplot(aes(x = MaxMean, y = p_value_adj, label = Genes, color = highlight, text = paste("Genes:", Genes, "<br>x:", MaxMean, "<br>y:", p_value))) + 
        geom_point(size=3) +
        scale_color_manual(values = c("Highlighted" = "red", "Normal" = "blue")) +
        labs(x=paste0("MaxMean"),
             y=paste0("p_value")) +
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
    
    output$heatmap_plot <- renderPlot({
      req(selected_variable_genes())
      variable_genes <- selected_variable_genes()
      df <- dataset$filtered_data()
      gene_data <- df %>% select(-Symbols)
      heatmap_df <- gene_data %>%
        inner_join(variable_genes, by = "Genes") %>%
        select(-p_value, -p_value_adj, -MaxMean)
      

      gene_expression_z <- t(scale(t(heatmap_df[,2:ncol(heatmap_df)])))
      gene_expression_z[is.na(gene_expression_z)] <- 0 
      set.seed(123)
      
      if (input$cluster_options == 1) {
        # k means clustering
        kmeans_result <- kmeans(gene_expression_z, centers=input$clustering_k , iter.max=20, nstart=200)
        rownames(gene_expression_z) <- paste(rownames(gene_expression_z), "Cluster", kmeans_result$cluster, sep=" ")
        
        p <- pheatmap(gene_expression_z,
                      cluster_rows = FALSE,  # Disable hierarchical clustering of rows
                      cluster_cols = FALSE,  # Disable hierarchical clustering of columns
                      # color = colorRampPalette(c("blue", "white", "red"))(50),
                      show_rownames = TRUE,
                      show_colnames = TRUE,
                      main = "K-means Clustering Gene Expression Matrix",
                      treeheight_row = 30, 
                      angle_col = 315, 
                      fontsize_col = 7,
                      width = 10, 
                      height = 12)
        print(p)
        return(p)
      } else {
        # hierarchical clustering
        pp <- pheatmap(gene_expression_z,
                      cluster_rows = TRUE,
                      clustering_distance_rows = "correlation",
                      cluster_cols = FALSE, 
                      fontsize_row = .10,
                      main = paste0("Correlation Distance Hierarchical Clustering Gene Expression Matrix"),
                      treeheight_row = 30, 
                      angle_col = 315, 
                      fontsize_col = 7,
                      width = 10, 
                      height = 12)    
        print(pp)
        return(pp)
      }
      # print(p)
      # return(p)
    })
    

  })
}

clusteringTabUI <- function(id) {
  fluidPage(
    titlePanel("Clustering of Gene Expression"),
    fluidRow(
      column(2, numericInput(NS(id, "x_cutoff"), "X Axis Cutoff", value = 0)),
      column(2, numericInput(NS(id, "y_cutoff"), "Y Axis Cutoff", value = 0)),
      column(2, selectInput(NS(id, "cluster_options"), "Clustering Methods",
                            choices = c("k means" = 1, "hierarchical" = 2), selected = 1)),
      column(2, numericInput(NS(id, "clustering_k"), "K Values", value = 5))
    ),
    fluidRow(
      column(2, actionButton(NS(id, "run_analysis"), "Run Analysis"))
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
