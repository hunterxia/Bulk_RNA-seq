library(tidyverse)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(plotly)
library(dplyr)
library(heatmaply)
library(logger)
library(shinybusy)

clusteringTabServer <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {

    # Reactive values to store data
    selected_variable_genes <- reactiveVal()
    samples_anova_results <- reactiveVal()
    clusters_download_data <- reactiveVal()
    hierarchical_distance_matrix <- reactiveVal()

    # Apply filters based on user inputs
    observeEvent(input$run_analysis, {
      isolate({
        log_info("Running analysis")
        # showModal(modalDialog("Running analysis, please wait...", footer = NULL))
        # Only recalculates when the above reactive values change
        tryCatch({
          result <- calculate_variable_genes()
          if (is.null(result)) {
            log_warn("Experimental group has less than 2 levels")
            # removeModal()
            showModal(modalDialog(
              title = "Error",
              "The factor 'EXPERIMENTAL_GROUP' has less than 2 levels. Please check your data.",
              footer = tagList(
                actionButton(NS(id, "close_modal"), "Close", class = "btn btn-primary")
              )))
          }
          else {
            log_info("ANOVA completed: {nrow(result)} genes processed")
            samples_anova_results(result)
          }
        }, error = function(e) {
          log_error("Analysis failed: {e$message}")
          stop(e)
        })
        # result <- calculate_variable_genes()
        # if (is.null(result)) {
        #   removeModal()
        #   showModal(modalDialog(
        #     title = "Error",
        #     "The factor 'EXPERIMENTAL_GROUP' has less than 2 levels. Please check your data.",
        #     footer = tagList(
        #       actionButton(NS(id, "close_modal"), "Close", class = "btn btn-primary")
        #     )))
        # }
        # else {
        #   samples_anova_results(result)
        # }
      })

    })

    # Close warning modal
    observeEvent(input$close_modal, {
      removeModal()
    })

    # calculate the ANOVA value
    calculate_variable_genes <- function() {

      req(dataset$filtered_data(), dataset$groups_data)

      df <- dataset$filtered_data()
      groups <- dataset$groups_data()

      # remove this column for simplicity 
      gene_data <- df %>% select(-Gene_Symbol)
      colnames(gene_data)[1] <- "Symbol"
      colnames(groups)[1] <- "Condition"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"


      gene_data_long <- gene_data %>%
        pivot_longer(
          cols = -Symbol,  # Exclude Gene_Symbol from the melting process
          names_to = "Condition",  # This will contain the sample/condition names
          values_to = "Expression"  # This will contain the expression values
        )

      # Join the dataframes to map Exp_Grp values to Conditions in gene_data_long
      gene_data_long <- gene_data_long %>%
        left_join(groups, by = "Condition")

      # Check if EXP_GROUP has at least 2 levels
      if (length(levels(factor(gene_data_long$EXPERIMENTAL_GROUP))) < 2) {
        return(NULL)
      }

      group_max <- gene_data_long %>%
        group_by(Symbol, EXPERIMENTAL_GROUP) %>%
        summarise(MaxExpression = max(Expression, na.rm = TRUE), .groups = 'drop')

      max_group_expression <- group_max %>%
        group_by(Symbol) %>%
        summarise(Max = max(MaxExpression))


      anova_results <- gene_data_long %>%
        group_by(Symbol) %>%
        do({
          model <- lm(Expression ~ EXPERIMENTAL_GROUP, data = .)
          anova_model <- anova(model)
          data.frame(p_value = anova_model$`Pr(>F)`[1])
        }) %>%
        ungroup()

      anova_results$p_value_adj <- p.adjust(anova_results$p_value, method = "BH")
      anova_results$p_value_log10 <- -log10(anova_results$p_value)
      anova_results %>%
        left_join(max_group_expression, by = "Symbol") %>%
        left_join(df %>% select(Symbol, Gene_Symbol), by = "Symbol")
    }

    # calculate the LFC value
    calculate_variable_genes_by_LFC <- function() {
      req(input$LFC_group_1, input$LFC_group_2)

      gene_data <- dataset$filtered_data()
      groups <- dataset$groups_data()
      selected_groups <- dataset$selected_groups()

      # remove this column for simplicity 
      #gene_data <- gene_data %>% select(-Symbol)
      colnames(gene_data)[1] <- "Symbol"
      colnames(gene_data)[2] <- "Gene_Symbol"
      colnames(groups)[1] <- "Condition"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"

      calculateGroupMeans <- function(normalized_counts, sample_info, selected_groups) {
        # Create an empty dataframe to store the results
        averaged_counts <- data.frame(row.names = rownames(normalized_counts))

        # Loop through each unique experimental group in the sample information
        for (Group in selected_groups) {
          # Find the samples belonging to this experimental group
          samples_in_grp <- sample_info$Condition[sample_info$EXPERIMENTAL_GROUP == Group]

          # Subset the normalized counts matrix to only these samples
          counts_subset <- normalized_counts[, colnames(normalized_counts) %in% samples_in_grp, drop = FALSE]

          # Calculate the mean across the columns (samples) for each gene
          group_means <- rowMeans(counts_subset, na.rm = TRUE)

          # Add the results to the dataframe, naming the column after the experimental group
          averaged_counts[[Group]] <- group_means
        }

        return(averaged_counts)
      }

      # apply group means function
      averaged_counts_result <- calculateGroupMeans(gene_data, groups, selected_groups)

      averaged_counts_result$max_mean <- apply(averaged_counts_result[, -1], 1, max)

      averaged_counts_result$max_mean_log2 <- log2(averaged_counts_result$max_mean + 1)

      averaged_counts_result <- averaged_counts_result %>%
        mutate(
          max_group = pmax(!!sym(input$LFC_group_1), !!sym(input$LFC_group_2)),
          min_group = pmin(!!sym(input$LFC_group_1), !!sym(input$LFC_group_2)),
          LFC = log2(max_group + 1) - log2(min_group + 1)  
        ) %>% 
        select(-max_group, -min_group)
      

      averaged_counts_result$Gene_Symbol <- gene_data$Gene_Symbol

      averaged_counts_result$Symbol <- gene_data$Symbol

      averaged_counts_result$LFC <- abs(averaged_counts_result$LFC)

      averaged_counts_result
    }

    # render variable genes plot
    # output$variable_genes_plot <- renderPlotly({
    #   on.exit(removeModal())
    #
    #   # render by ANOVA
    #   if (input$y_axis_var == 1) {
    #     req(samples_anova_results())
    #     result_data <- samples_anova_results()
    #     result_data$Max_log2 <- log2(result_data$Max + 1)
    #
    #     # update selected variable genes
    #     variable_genes <- dplyr::filter(result_data, p_value <= input$y_cutoff, Max_log2 >= input$x_cutoff)
    #     selected_variable_genes(variable_genes)
    #
    #     # highlight selected variable genes in the plot
    #     result_data$highlight <- ifelse(result_data$p_value <= input$y_cutoff & result_data$Max_log2 >= input$x_cutoff, "Selected", "Unselected")
    #
    #     p <- result_data %>%
    #       ggplot(aes(x = Max_log2, y = p_value_log10, label = Gene_Symbol, color = highlight,
    #                  text = paste("Genes:", Gene_Symbol, "<br>Max Group Expression:", Max, "<br>P Value:", p_value, "<br>P Adjusted Value:", p_value_adj, "<br>-log10(P-value):", p_value_log10))) +
    #       geom_point(size = 3) +
    #       scale_color_manual(values = c("Selected" = "blue", "Unselected" = "grey")) +
    #       labs(x = paste0("log2(Max) + 1"),
    #            y = paste0("-log10(P-value)")) +
    #       theme_minimal() +
    #       theme(legend.position = "right") +
    #       theme_linedraw(base_size = 16) +
    #       theme(panel.grid.major = element_blank(),
    #             panel.grid.minor = element_blank(),
    #             legend.title = element_blank(),
    #             legend.text = element_text(size = 10),
    #             legend.position = "right")
    #   } else {
    #     # render by LFC
    #     if (input$LFC_group_1 == input$LFC_group_2) {
    #       return()
    #     }
    #
    #     result_data <- calculate_variable_genes_by_LFC()
    #     # update selected variable genes
    #     variable_genes <- dplyr::filter(result_data, LFC >= input$y_cutoff, max_mean_log2 >= input$x_cutoff)
    #     selected_variable_genes(variable_genes)
    #
    #     # highlight selected variable genes in the plot
    #     result_data$highlight <- ifelse(result_data$LFC >= input$y_cutoff & result_data$max_mean_log2 >= input$x_cutoff, "Selected", "Unselected")
    #     p <- result_data %>%
    #       ggplot(aes(x = max_mean_log2, y = LFC, color = highlight, text = paste("Genes:", Gene_Symbol, "<br>log2(Max-mean) + 1:", max_mean_log2, "<br>LFC:", LFC))) +
    #       geom_point(size = 3) +
    #       scale_color_manual(values = c("Selected" = "blue", "Unselected" = "grey")) +
    #       labs(x = paste0("log2(Max-mean) + 1"),
    #            y = paste0("LFC")) +
    #       theme_minimal() +
    #       theme(legend.position = "right") +
    #       theme_linedraw(base_size = 16) +
    #       theme(panel.grid.major = element_blank(),
    #             panel.grid.minor = element_blank(),
    #             legend.title = element_blank(),
    #             legend.text = element_text(size = 10),
    #             legend.position = "right")
    #
    #   }
    #   ggplotly(p, tooltip = "text")
    # })

    output$variable_genes_plot <- renderPlotly({
      # on.exit(removeModal())

      # render by ANOVA
      if (input$y_axis_var == 1) {
        req(samples_anova_results())
        result_data <- samples_anova_results()
        result_data$Max_log2 <- log2(result_data$Max + 1)

        variable_genes <- dplyr::filter(result_data, p_value <= input$y_cutoff, Max_log2 >= input$x_cutoff)
        selected_variable_genes(variable_genes)

        result_data$highlight <- ifelse(
          result_data$p_value <= input$y_cutoff & result_data$Max_log2 >= input$x_cutoff,
          "Selected", "Unselected"
        )

        p <- plot_ly(
          data = result_data,
          x = ~Max_log2,
          y = ~p_value_log10,
          type = 'scatter',
          mode = 'markers',
          color = ~highlight,
          colors = c("Selected" = "blue", "Unselected" = "grey"),
          text = ~paste(
            "Genes:", Gene_Symbol,
            "<br>Max Group Expression:", Max,
            "<br>P Value:", p_value,
            "<br>P Adjusted Value:", p_value_adj,
            "<br>-log10(P-value):", p_value_log10
          ),
          hoverinfo = "text",
          marker = list(size = 10)
        ) %>%
          layout(
            title = "Variable Genes Plot (ANOVA)",
            xaxis = list(title = "log2(Max + 1)"),
            yaxis = list(title = "-log10(P-value)"),
            legend = list(orientation = "v", x = 1, y = 1)
          )

      } else {
        # render by LFC
        if (input$LFC_group_1 == input$LFC_group_2) {
          return()
        }

        result_data <- calculate_variable_genes_by_LFC()
        variable_genes <- dplyr::filter(result_data, LFC >= input$y_cutoff, max_mean_log2 >= input$x_cutoff)
        selected_variable_genes(variable_genes)

        result_data$highlight <- ifelse(
          result_data$LFC >= input$y_cutoff & result_data$max_mean_log2 >= input$x_cutoff,
          "Selected", "Unselected"
        )

        p <- plot_ly(
          data = result_data,
          x = ~max_mean_log2,
          y = ~LFC,
          type = 'scatter',
          mode = 'markers',
          color = ~highlight,
          colors = c("Selected" = "blue", "Unselected" = "grey"),
          text = ~paste(
            "Genes:", Gene_Symbol,
            "<br>log2(Max-mean + 1):", max_mean_log2,
            "<br>LFC:", LFC
          ),
          hoverinfo = "text",
          marker = list(size = 10)
        ) %>%
          layout(
            title = "Variable Genes Plot (LFC)",
            xaxis = list(title = "log2(Max+1)"), 
            yaxis = list(title = "LFC"),
            legend = list(orientation = "v", x = 1, y = 1)
          )
      }

      p
    })

    output$number_of_genes <- renderText(paste0("Number of Selected Genes: ", nrow(selected_variable_genes())))

    # render heatmap using selected variable genes
    create_cluster_plot <- function(data, cluster_options, max_clusters) {
      data_sample_expr <- data[, -c(1, 2)]

      min_max_normalize <- function(x) {
        (x - min(x)) / (max(x) - min(x))
      }

      gene_expression_z <- t(apply(data_sample_expr, 1, min_max_normalize))
      gene_expression_z[is.na(gene_expression_z)] <- 0
      rownames(gene_expression_z) <- data$Symbol

      if (cluster_options == 1) {
        if (nrow(data_sample_expr) < max_clusters) {
          showNotification("Error: Number of data points is less than the number of clusters chosen.", type = "error")
          return(NULL)
        }
        set.seed(40)
        clusters <- kmeans(gene_expression_z, centers = max_clusters, nstart = 25)

        # Create annotation dataframe with Gene_Symbol and Clusters
        annotation_df <- data.frame(Gene_Symbol = data$Gene_Symbol, Clusters = clusters[["cluster"]])
        rownames(annotation_df) <- data$Symbol
        
        # Create a more comprehensive dataframe for download that includes expression values
        download_df <- data.frame(
          Symbol = data$Symbol,
          Gene_Symbol = data$Gene_Symbol,
          Cluster = clusters[["cluster"]]
        )
        
        # Add all expression values from the original data
        download_df <- cbind(download_df, data_sample_expr)
        
        # Store the comprehensive dataframe for download
        clusters_download_data(download_df)
        
        # Use the simpler annotation dataframe for the heatmap display
        annotation_df <- annotation_df %>%
          arrange(Clusters) %>%
          select(-Gene_Symbol)
        order <- order(clusters$cluster)
        cor_matrix_ordered <- gene_expression_z[order,]

        # Basic interactive heatmap with simplified hover functionality
        gene_labels <- data$Gene_Symbol[order]
        sample_labels <- colnames(data_sample_expr)
        
        plot_data <- as.data.frame(cor_matrix_ordered)
        plot_data$Gene <- gene_labels
        plot_data_long <- pivot_longer(plot_data, cols = -Gene, names_to = "Sample", values_to = "Normalized_Value")
        # Ensure no duplicates
        plot_data_long <- distinct(plot_data_long, Gene, Sample, .keep_all = TRUE)
        
        # Add cluster information
        plot_data_long$Cluster <- as.factor(clusters$cluster[order][match(plot_data_long$Gene, gene_labels)])
        
        # Add the original expression values for hover
        original_values <- as.data.frame(as.matrix(data_sample_expr)[order,])
        original_values$Gene <- gene_labels
        original_values_long <- pivot_longer(original_values, cols = -Gene, names_to = "Sample", values_to = "Expression_Value")
        # Ensure no duplicates
        original_values_long <- distinct(original_values_long, Gene, Sample, .keep_all = TRUE)
        
        # Combine normalized and original values
        plot_data_final <- left_join(plot_data_long, original_values_long, by = c("Gene", "Sample"), relationship = "one-to-one")
        
        # Create the plotly heatmap
        p <- plot_ly(
          data = plot_data_final,
          x = ~Sample,
          y = ~Gene,
          z = ~Normalized_Value,
          type = "heatmap",
          colors = colorRamp(c("navy", "white", "firebrick3")),
          text = ~paste(
            "Gene:", Gene,
            "<br>Sample:", Sample,
            "<br>Expression Value:", round(Expression_Value, 2),
            "<br>Normalized Value:", round(Normalized_Value, 4),
            "<br>Cluster:", Cluster
          ),
          hoverinfo = "text"
        ) %>%
        layout(
          title = "K-means Clustered Matrix",
          xaxis = list(title = "", tickangle = 45),
          yaxis = list(title = "")
        )

        return(p)
      }

      # using correlation to calculate distance matrix
      row_cor_matrix <- cor(t(gene_expression_z))
      row_distance_matrix <- as.dist(1 - row_cor_matrix)
      
      # Create a comprehensive dataframe for hierarchical clustering download
      # This includes the distance matrix and expression values
      hier_download_df <- as.data.frame(1 - row_cor_matrix)
      hier_download_df$Symbol <- data$Symbol
      hier_download_df$Gene_Symbol <- data$Gene_Symbol
      
      # Add expression values
      expression_data <- as.data.frame(gene_expression_z)
      hier_download_df <- cbind(hier_download_df[c("Symbol", "Gene_Symbol")], expression_data, hier_download_df[!(names(hier_download_df) %in% c("Symbol", "Gene_Symbol"))])
      
      hierarchical_distance_matrix(hier_download_df)
      
      # Create a simplified interactive hierarchical heatmap
      hclust_rows <- hclust(row_distance_matrix, method = "complete")
      
      # Reorder data by hierarchical clustering
      hc_order <- hclust_rows$order
      gene_labels <- data$Gene_Symbol[hc_order]
      
      # Prepare the long-format data for plotly
      plot_data <- as.data.frame(gene_expression_z[hc_order,])
      plot_data$Gene <- gene_labels
      plot_data_long <- pivot_longer(plot_data, cols = -Gene, names_to = "Sample", values_to = "Normalized_Value")
      # Ensure no duplicates
      plot_data_long <- distinct(plot_data_long, Gene, Sample, .keep_all = TRUE)
      
      # Add the original expression values for hover
      original_values <- as.data.frame(as.matrix(data_sample_expr)[hc_order,])
      original_values$Gene <- gene_labels
      original_values_long <- pivot_longer(original_values, cols = -Gene, names_to = "Sample", values_to = "Expression_Value")
      # Ensure no duplicates
      original_values_long <- distinct(original_values_long, Gene, Sample, .keep_all = TRUE)
      
      # Combine normalized and original values
      plot_data_final <- left_join(plot_data_long, original_values_long, by = c("Gene", "Sample"), relationship = "one-to-one")
      
      # Create the plotly heatmap
      p <- plot_ly(
        data = plot_data_final,
        x = ~Sample,
        y = ~Gene,
        z = ~Normalized_Value,
        type = "heatmap",
        colors = colorRamp(c("navy", "white", "firebrick3")),
        text = ~paste(
          "Gene:", Gene,
          "<br>Sample:", Sample,
          "<br>Expression Value:", round(Expression_Value, 2),
          "<br>Normalized Value:", round(Normalized_Value, 4)
        ),
        hoverinfo = "text"
      ) %>%
      layout(
        title = "Hierarchically Clustered Matrix",
        xaxis = list(title = "", tickangle = 45),
        yaxis = list(title = "")
      )

      return(p)
    }

    output$heatmap_plot <- renderPlotly({
      req(dataset$filtered_data(), selected_variable_genes())
      variable_genes <- selected_variable_genes()
      gene_data <- dataset$filtered_data()
      groups <- dataset$groups_data()
      variable_genes <- variable_genes %>% select(Symbol)

      # show by group
      if (input$clustering_grouped) {
        gene_data_long <- pivot_longer(gene_data, cols = -names(gene_data)[1:2], names_to = "Sample", values_to = "Values")

        gene_data_long_grouped <- gene_data_long %>%
          left_join(groups, by = "Sample")

        gene_data_long_grouped_sum <- gene_data_long_grouped %>%
          group_by(Symbol, Gene_Symbol, Group) %>%
          summarize(Values = sum(Values), .groups = "drop")

        gene_data <- pivot_wider(gene_data_long_grouped_sum, names_from = Group, values_from = Values, values_fill = list(value = 0))
      }

      # by default, show by each sample
      heatmap_df <- gene_data %>%
        inner_join(variable_genes, by = "Symbol")
      create_cluster_plot(heatmap_df, input$cluster_options, input$clustering_k)
    })

    # download k-means clustering data
    output$download_clusters <- downloadHandler(
      filename = function() {
        paste("assigned-clusters-with-expression-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(clusters_download_data(), file, row.names = FALSE)
      }
    )

    # download hierarchical clustering data
    output$download_hierarchical <- downloadHandler(
      filename = function() {
        paste("hierarchical-clustering-with-expression-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(hierarchical_distance_matrix(), file, row.names = FALSE)
      }
    )

    # Reading and storing expression data
    observeEvent(input$y_axis_var, {
      groups <- dataset$selected_groups()
      options <- c()
      for (group in groups) {
        options <- c(options, group)
      }
      
      output$y_axis_operations <- renderUI({
        if(input$y_axis_var == 2) {
          fluidRow(
            column(6,
                   selectInput(session$ns("LFC_group_1"), "Group 1", options, selected = options[1])
            ),
            column(6,
                   selectInput(session$ns("LFC_group_2"), "Group 2", options, selected = options[2])
            )
          )
        } else {
          # Return empty UI for ANOVA since we don't need additional controls
          NULL
        }
      })
      
      output$y_cutoff <- renderUI({
        if(input$y_axis_var == 1) {
          tagList(
            numericInput(NS(id, "y_cutoff"), "P-value Cutoff", value = 0.05, min = 0, max = 1, step = 0.1),
            actionButton(NS(id, "run_analysis"), "Run Analysis")
          )
        } else {
          numericInput(NS(id, "y_cutoff"), "Fold Change Cutoff", value = 1, step = 0.1,min = 0)
        }
      })
    })
    
    observeEvent(input$cluster_options, {
      if (input$cluster_options == 1) {
        output$clustering_operations <- renderUI({
          fluidRow(
            div(class = "bottom-centered",
                column(5, numericInput(NS(id, "clustering_k"), "K Values", value = 5, step = 1)),
                column(6, downloadButton(NS(id, "download_clusters"), "Download clusters with expression"))
            )
          )
        })

      } else {
        output$clustering_operations <- renderUI({
          fluidRow(
            div(class = "bottom-centered",
                column(6, downloadButton(NS(id, "download_hierarchical"), "Download hierarchical data with expression"))
            )
          )
        })
      }
    })

    observeEvent(input$LFC_group_1, {
      if (input$LFC_group_1 == input$LFC_group_2) {
        showNotification("Please select two different gourps for LFC calculation.", type = "error")
      }
    })

    observeEvent(input$LFC_group_2, {
      if (input$LFC_group_1 == input$LFC_group_2) {
        showNotification("Please select two different gourps for LFC calculation.", type = "error")
      }
    })

    # Return reactive values
    return(list(
      selected_variable_genes = selected_variable_genes
    ))


  })
}

clusteringTabUI <- function(id) {
  fluidPage(
    add_busy_spinner(spin = "fading-circle", color = "#000000", position = "top-right"),
    tags$head(
      tags$style(HTML("
        .bottom-centered {
          display: flex;
          flex-direction: column;
          align-items: flex-start;
          gap: 10px;
        }
        .shiny-notification {
          position: fixed;
          top: 50% !important;
          left: 50% !important;
          transform: translate(-50%, -50%);
          font-size: 16px;
          color: white;
          background-color: red;
          padding: 10px;
          border-radius: 5px;
        }
        .js-plotly-plot {
          margin: 0 auto;
        }
        .plotly-heatmap {
          width: 100% !important;
          height: 100% !important;
        }
        .modebar {
          background-color: rgba(255, 255, 255, 0.8) !important;
        }
      "))
    ),
    titlePanel("Clustering of Gene Expression"),
    
    # Step 1: Variable Gene Selection
    fluidRow(
      column(12, 
             h3("Select Variable Genes"),
             hr()
      )
    ),
    fluidRow(
      column(4,
             div(class = "bottom-centered",
                 selectInput(NS(id, "y_axis_var"), "Y Axis Variability",
                             c("ANOVA" = 1, "LFC" = 2), selected = 1),
                 uiOutput(NS(id, "y_axis_operations")),
                 numericInput(NS(id, "x_cutoff"), "log2(Max+1) Cutoff", value = 1, min = 0, step = 0.1),
                 uiOutput(NS(id, "y_cutoff"))
             )
      ),
      column(8,
             textOutput(NS(id, "number_of_genes")),
             plotlyOutput(NS(id, "variable_genes_plot"), width = "100%", height = "50vh")
      )
    ),
    
    # Step 2: Clustering
    fluidRow(
      column(12, 
             h3("Cluster Selected Genes"),
             hr()
      )
    ),
    fluidRow(
      column(4,
             div(class = "bottom-centered",
                 materialSwitch(inputId = NS(id, "clustering_grouped"), 
                                label = "Cluster by group:", 
                                value = FALSE, 
                                status = "primary"),
                 selectInput(NS(id, "cluster_options"), "Clustering Method",
                             choices = c("k means" = 1, "hierarchical" = 2), selected = 1),
                 uiOutput(NS(id, "clustering_operations"))
             )
      ),
      column(8,
             plotlyOutput(NS(id, "heatmap_plot"), width = "100%", height = "60vh")
      )
    )
  )
}
