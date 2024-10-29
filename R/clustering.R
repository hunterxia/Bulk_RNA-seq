library(tidyverse)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(plotly)
library(shinybusy)
library(dplyr)
#library(rlang)

clusteringTabServer <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    # Reactive values to store data
    selected_variable_genes <- reactiveVal()
    samples_anova_results <- reactiveVal()
    clusters_download_data <- reactiveVal()
    
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
      gene_data <- df %>% select(-Symbol)
      colnames(gene_data)[1] <- "Gene_Symbol"
      colnames(groups)[1] <- "Condition"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"
      
      
      gene_data_long <- gene_data %>%
        pivot_longer(
          cols = -Gene_Symbol,  # Exclude Gene_Symbol from the melting process
          names_to = "Condition",  # This will contain the sample/condition names
          values_to = "Expression"  # This will contain the expression values
        )
      
      # Join the dataframes to map Exp_Grp values to Conditions in gene_data_long
      gene_data_long <- gene_data_long %>%
        left_join(groups, by = "Condition")
      
      group_max <- gene_data_long %>%
        group_by(Gene_Symbol, EXPERIMENTAL_GROUP) %>%
        summarise(MaxExpression = max(Expression, na.rm = TRUE), .groups = 'drop')
      
      max_group_expression <- group_max %>%
        group_by(Gene_Symbol) %>%
        summarise(Max = max(MaxExpression))
      
      
      anova_results <- gene_data_long %>%
        group_by(Gene_Symbol) %>%
        do({
          model <- lm(Expression ~ EXPERIMENTAL_GROUP, data = .)
          anova_model <- anova(model)
          data.frame(p_value = anova_model$`Pr(>F)`[1])
        }) %>%
        ungroup()
      
      anova_results$p_value_adj <- p.adjust(anova_results$p_value, method = "BH")
      anova_results$p_value_log10 <- -log10(anova_results$p_value)
      anova_results %>%
        left_join(max_group_expression, by = "Gene_Symbol")
      
    }
    
    calculate_variable_genes_by_LFC <- function() {
      req(input$LFC_group_1, input$LFC_group_2)
      
      gene_data <- dataset$filtered_data()
      groups <- dataset$groups_data()
      selected_groups <- dataset$selected_groups()

      # remove this column for simplicity 
      gene_data <- gene_data %>% select(-Symbol)
      colnames(gene_data)[1] <- "Gene_Symbol"
      colnames(groups)[1] <- "Condition"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"
      
      calculateGroupMeans <- function(normalized_counts, sample_info, selected_groups) {
        # Create an empty dataframe to store the results
        averaged_counts <- data.frame(row.names = rownames(normalized_counts))
        
        # Loop through each unique experimental group in the sample information
        for(Group in selected_groups) {
          # Find the samples belonging to this experimental group
          samples_in_grp <- sample_info$Condition[sample_info$EXPERIMENTAL_GROUP == Group]
          
          # Subset the normalized counts matrix to only these samples
          counts_subset <- normalized_counts[, colnames(normalized_counts) %in% samples_in_grp, drop = FALSE]
          
          # Calculate the mean across the columns (samples) for each gene
          group_means <- rowMeans(counts_subset, na.rm = TRUE)
          
          # Add the results to the dataframe, naming the column after the experimental group
          averaged_counts[[exp_grp]] <- group_means
        }
        
        return(averaged_counts)
      }
      
      # apply group means function
      averaged_counts_result <- calculateGroupMeans(gene_data, groups, selected_groups)

      averaged_counts_result$max_mean <- apply(averaged_counts_result[, -1], 1, max)
      
      averaged_counts_result$max_mean_log2 <- log2(averaged_counts_result$max_mean + 1)
      
      averaged_counts_result <- averaged_counts_result %>%
        mutate(LFC = log2(!!sym(input$LFC_group_1) + 1e-4) - log2(!!sym(input$LFC_group_2) + 1e-4))
      
      averaged_counts_result$Gene_Symbol <- gene_data$Gene_Symbol
      
      averaged_counts_result$LFC <- abs(averaged_counts_result$LFC)
      
      averaged_counts_result
    }
    
    output$variable_genes_plot <- renderPlotly({
      if (input$y_axis_var == 1) {
        req(samples_anova_results())
        result_data <- samples_anova_results()
  
        result_data$Max_log2 <- log2(result_data$Max + 1)
  
        # update selected variable genes
        variable_genes <- dplyr::filter(result_data, p_value_log10 <= input$y_cutoff, Max_log2 >= input$x_cutoff)
        selected_variable_genes(variable_genes)
        
        # highlight selected variable genes in the plot
        result_data$highlight <- ifelse(result_data$p_value_log10 <= input$y_cutoff & result_data$Max_log2 >= input$x_cutoff, "Selected", "Unselected")
        
        p <- result_data %>%
          ggplot(aes(x = Max_log2, y = p_value_log10, label = Gene_Symbol, color = highlight,
                     text = paste("Genes:", Gene_Symbol, "<br>Max Group Expression:", Max, "<br>P Value:", p_value, "<br>P Adjusted Value:", p_value_adj, "<br>P Value(-log10):", p_value_log10))) + 
          geom_point(size=3) +
          scale_color_manual(values = c("Selected" = "blue", "Unselected" = "grey")) +
          labs(x=paste0("Max Group Expression (log2 + 1)"),
               y=paste0("P Values (-log10)")) +
          theme_minimal() + 
          theme(legend.position = "right") +
          theme_linedraw(base_size=16) +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                legend.title = element_blank(),
                legend.text = element_text(size = 10),
                legend.position = "right")
      } else {
        result_data <- calculate_variable_genes_by_LFC()
        
        # update selected variable genes
        variable_genes <- dplyr::filter(result_data, LFC >= input$y_cutoff, max_mean_log2 >= input$x_cutoff)
        selected_variable_genes(variable_genes)
        
        # highlight selected variable genes in the plot
        result_data$highlight <- ifelse(result_data$LFC >= input$y_cutoff & result_data$max_mean_log2 >= input$x_cutoff, "Selected", "Unselected")
        p <- result_data %>%
          ggplot(aes(x = max_mean_log2, y = LFC, color = highlight, text = paste("Genes:", Gene_Symbol, "<br>max mean by log2+1:", max_mean_log2, "<br>LFC:", LFC))) + 
          geom_point(size=3) +
          scale_color_manual(values = c("Selected" = "blue", "Unselected" = "grey")) +
          labs(x=paste0("Max Mean (log2 + 1)"),
               y=paste0("LFC")) +
          theme_minimal() + 
          theme(legend.position = "right") +
          theme_linedraw(base_size=16) +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                legend.title = element_blank(),
                legend.text = element_text(size = 10),
                legend.position = "right")

      }
      ggplotly(p, tooltip = "text")
    })
    
    output$number_of_genes <- renderText(paste0("Number of Selescted Genes: ", nrow(selected_variable_genes())))
    
    create_cluster_plot <- function(data, cluster_options, max_clusters) {
      data_sample_expr <- data[,-c(1, 2)]
      #gene_expression_z <- scale(data_sample_expr)
      
      min_max_normalize <- function(x) {
        (x - min(x)) / (max(x) - min(x))
      }

      gene_expression_z <- t(apply(data_sample_expr, 1, min_max_normalize))
      gene_expression_z[is.na(gene_expression_z)] <- 0 
      rownames(gene_expression_z) <- data$Symbol

      if (cluster_options == 1) {
        set.seed(40)
        clusters <- kmeans(gene_expression_z, centers=max_clusters, nstart=25)

        annotation_df <- data.frame(Gene_Symbol = data$Gene_Symbol, Clusters = clusters[["cluster"]])
        rownames(annotation_df) <- data$Symbol
        annotation_df <- annotation_df %>% arrange(Clusters) %>% select(-Gene_Symbol)
        clusters_download_data(annotation_df)
        order <- order(clusters$cluster)
        cor_matrix_ordered <- gene_expression_z[order, ]
        p <- pheatmap::pheatmap(cor_matrix_ordered, annotation_row = annotation_df,
                                                               cluster_rows = FALSE,
                                                               cluster_cols = FALSE,
                                                               fontsize_row = 10,
                                                               main = "K-means Clustered Matrix",
                                                               treeheight_row = 30,
                                                               labels_row = "",
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
      # gene_data <- df %>% select(-Symbol)
      gene_data <- df
      variable_genes <- variable_genes %>% select(Gene_Symbol)

      heatmap_df <- gene_data %>%
        inner_join(variable_genes, by = "Gene_Symbol")
      p <- create_cluster_plot(heatmap_df, input$cluster_options, input$clustering_k)
      return(p)
    })
    
    # download correlation data
    output$download_clusters <- downloadHandler(
      filename = function() {
        paste("assigned-clusters-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(clusters_download_data(), file)
      }
    )
    
    # Reading and storing expression data
    observeEvent(input$y_axis_var, {
      if (input$y_axis_var == 2) {
        groups <- dataset$selected_groups()
        options <- c()
        for (group in groups) {
          options <- c(options, group)
        }
        output$y_axis_operations <- renderUI({
           fluidRow(
             column(6, 
                    selectInput(session$ns("LFC_group_1"), "Group 1", options, selected=options[1])
             ),
             column(6, 
                    selectInput(session$ns("LFC_group_2"), "Group 2", options, selected=options[2])
             ),
           )
        })
        
        output$y_cutoff <- renderUI({
          numericInput(NS(id, "y_cutoff"), "LFC Cutoff", value = 1)
        })
      } else {
        output$y_axis_operations <- renderUI({
          actionButton(NS(id, "run_analysis"), "Run Analysis")
        })
        output$y_cutoff <- renderUI({
          numericInput(NS(id, "y_cutoff"), "P Values(-log10) Cutoff", value = 1)
        })
      }
    })
    
    
  })
}

clusteringTabUI <- function(id) {
  fluidPage(
    tags$head(
      tags$style(HTML("
      .bottom-centered {
        display: flex;
        align-items: center;
      }
    "))
    ),
    titlePanel("Clustering of Gene Expression"),
    fluidRow(class = "bottom-centered",
      column(2, 
             selectInput(NS(id, "y_axis_var"), "Y Axis Variability",
                         c("ANOVA" = 1, "LFC" = 2), selected=1)
      ),
      column(4,  uiOutput(NS(id, ("y_axis_operations"))))
    ),
    fluidRow(
      div(class = "bottom-centered",
        column(2, numericInput(NS(id, "x_cutoff"), "Max Group Expression(log2+1) Cutoff", value = 0)),
        column(2, uiOutput(NS(id, ("y_cutoff")))),
        column(2, selectInput(NS(id, "cluster_options"), "Clustering Methods",
                              choices = c("k means" = 1, "hierarchical" = 2), selected = 1)),
        column(2, numericInput(NS(id, "clustering_k"), "K Values", value = 5)),
        column(2, downloadButton(NS(id, "download_clusters"), "Download assigned clusters"))
      )
    ),
    fluidRow(
      column(6,
             textOutput(NS(id, "number_of_genes")),
             add_busy_spinner(spin = "fading-circle", color = "#000000"),
             plotlyOutput(NS(id, "variable_genes_plot"), width = "100%", height = "700px")),
      column(6, plotOutput(NS(id, "heatmap_plot"), width = "100%", height = "700px"))
    ),
  )
  
}
