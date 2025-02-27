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

        annotation_df <- data.frame(Gene_Symbol = data$Gene_Symbol, Clusters = clusters[["cluster"]])
        rownames(annotation_df) <- data$Symbol
        clusters_download_data(annotation_df)
        annotation_df <- annotation_df %>%
          arrange(Clusters) %>%
          select(-Gene_Symbol)
        order <- order(clusters$cluster)
        cor_matrix_ordered <- gene_expression_z[order,]

        # built-in heatmap has high resolution, but does not have move over function
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

      # using correlation to calculate distance matrix
      row_cor_matrix <- cor(t(gene_expression_z))
      row_distance_matrix <- as.dist(1 - row_cor_matrix)
      hierarchical_distance_matrix(1 - row_cor_matrix)
      p <- pheatmap::pheatmap(gene_expression_z,
                              cluster_rows = TRUE,
                              cluster_cols = FALSE,
                              clustering_distance_rows = row_distance_matrix,
                              fontsize_row = 10,
                              main = "Hierarchically Clustered Matrix",
                              treeheight_row = 30,
                              angle_col = 315,
                              show_rownames = FALSE,
                              fontsize_col = 10,
                              width = 10,
                              height = 10)

      return(p)
    }

    output$heatmap_plot <- renderPlot({
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
      p <- create_cluster_plot(heatmap_df, input$cluster_options, input$clustering_k)
      return(p)
    })

    # download k-means clustering data
    output$download_clusters <- downloadHandler(
      filename = function() {
        paste("assigned-clusters-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(clusters_download_data(), file)
      }
    )

    # download hierarchical clustering data
    output$download_hierarchical <- downloadHandler(
      filename = function() {
        paste("hierarchical-clustering-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(hierarchical_distance_matrix(), file)
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
                column(5, numericInput(NS(id, "clustering_k"), "K Values", value = 5, step = 0.1)),
                # column(3, materialSwitch(inputId = NS(id, "clustering_grouped"), label = "Show by group: ", value = FALSE, status = "primary")),
                column(4, downloadButton(NS(id, "download_clusters"), "Download assigned clusters"))
            )
          )
        })

      } else {
        output$clustering_operations <- renderUI({
          fluidRow(
            div(class = "bottom-centered",
                # column(3,materialSwitch(inputId = NS(id, "clustering_grouped"), label = "Show by group: ", value = FALSE, status = "primary")),
                column(6, downloadButton(NS(id, "download_hierarchical"), "Download distance matrix"))
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
             plotOutput(NS(id, "heatmap_plot"), width = "100%", height = "50vh")
      )
    )
  )
}
