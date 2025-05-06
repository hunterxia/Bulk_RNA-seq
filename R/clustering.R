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
    variability_plot_data_cache <- reactiveVal(NULL) # Cache for variability plot base data

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

    # Calculate gene variability metric (Max-Min expression)
    calculate_gene_variability_metric <- function() {
      req(dataset$filtered_data())

      gene_data_full <- dataset$filtered_data() # This contains Symbol, Gene_Symbol, and sample columns

      # Ensure 'Symbol' and 'Gene_Symbol' are present
      if (!("Symbol" %in% names(gene_data_full)) || !("Gene_Symbol" %in% names(gene_data_full))) {
        log_error("Symbol or Gene_Symbol column missing in filtered_data")
        return(NULL)
      }
      
      # Select only sample columns for Min/Max calculation
      sample_columns <- setdiff(names(gene_data_full), c("Symbol", "Gene_Symbol"))
      
      if (length(sample_columns) == 0) {
        log_error("No sample columns found in filtered_data for variability calculation.")
        return(NULL)
      }
      
      # Calculate Min and Max expression for each gene across selected samples
      gene_expression_stats <- gene_data_full %>%
        rowwise() %>%
        mutate(
          Min_Expression = min(c_across(all_of(sample_columns)), na.rm = TRUE),
          Max_Expression = max(c_across(all_of(sample_columns)), na.rm = TRUE)
        ) %>%
        ungroup() %>%
        select(Symbol, Gene_Symbol, Min_Expression, Max_Expression)

      # Calculate Variability and overall average expression (for x-axis)
      variability_results <- gene_expression_stats %>%
        mutate(
          Variability = log2(Max_Expression + 1) - log2(Min_Expression + 1),
          Avg_Expression = rowMeans(gene_data_full[, sample_columns, drop = FALSE], na.rm = TRUE) # Calculate mean of original values
        ) %>%
        mutate(
          log2_Avg_Expression_plus_1 = log2(Avg_Expression + 1) # For the X-axis, consistent with previous "max_mean_log2"
        )

      # Handle cases where Variability might be NaN or Inf (e.g., if Min_Expression is -1 due to no expression before +1)
      variability_results <- variability_results %>%
        mutate(
          Variability = ifelse(is.finite(Variability), Variability, 0)
        )

      return(variability_results %>% select(Symbol, Gene_Symbol, log2_Avg_Expression_plus_1, Variability))
    }

    # Observer to update cached variability data when mode changes or main data changes
    observe({
      req(input$y_axis_var == 2, dataset$filtered_data()) # Only for Expression Variability mode
      log_info("Recalculating base data for variability plot due to data/mode change.")
      variability_plot_data_cache(calculate_gene_variability_metric())
    })

    output$variable_genes_plot <- renderPlotly({
      if (input$y_axis_var == 1) { # ANOVA
        req(samples_anova_results())
        result_data <- samples_anova_results()
        result_data$Max_log2 <- log2(result_data$Max + 1)

        current_x_cutoff_anova <- req(input$x_cutoff)
        current_y_cutoff_anova <- req(input$y_cutoff)

        variable_genes_anova <- dplyr::filter(result_data, p_value <= current_y_cutoff_anova, Max_log2 >= current_x_cutoff_anova)
        selected_variable_genes(variable_genes_anova)

        result_data$highlight <- ifelse(
          result_data$p_value <= current_y_cutoff_anova & result_data$Max_log2 >= current_x_cutoff_anova,
          "Selected", "Unselected"
        )

        p_anova <- plot_ly(
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
        ) %>% layout(
            title = "Variable Genes Plot (ANOVA)",
            xaxis = list(title = "log2(Max Expression + 1)"), # Clarified X-axis for ANOVA
            yaxis = list(title = "-log10(P-value)"),
            legend = list(orientation = "v", x = 1, y = 1, title = list(text = "Status"))
          )
        return(p_anova)

      } else { # Render by new Variability Metric (initial render)
        req(variability_plot_data_cache())
        plot_data <- variability_plot_data_cache()
        
        current_x_cutoff_var <- req(input$x_cutoff)
        current_y_cutoff_var <- req(input$y_cutoff)

        # Determine highlight status and direct colors for initial plot
        plot_data$highlight_color <- ifelse(
          plot_data$Variability >= current_y_cutoff_var & plot_data$log2_Avg_Expression_plus_1 >= current_x_cutoff_var,
          "blue", 
          "grey"  
        )
        
        current_selected_var_genes <- dplyr::filter(plot_data, Variability >= current_y_cutoff_var, log2_Avg_Expression_plus_1 >= current_x_cutoff_var)
        selected_variable_genes(current_selected_var_genes)

        p_var <- plot_ly(
          data = plot_data, 
          x = ~log2_Avg_Expression_plus_1, 
          y = ~Variability, 
          type = 'scatter', 
          mode = 'markers',
          marker = list(size = 10, color = ~highlight_color), # Use direct colors
          text = ~paste(
            "Gene:", Gene_Symbol,
            "<br>log2(Avg Expression + 1):", round(log2_Avg_Expression_plus_1, 2),
            "<br>Variability (log2(Max+1)-log2(Min+1)):", round(Variability, 2)
          ),
          hoverinfo = "text"
        ) %>% layout(
            title = "Variable Genes Plot (Expression Variability)",
            xaxis = list(title = "log2(Average Expression + 1)"), 
            yaxis = list(title = "Variability [log2(Max+1) - log2(Min+1)]"),
            legend = list(showlegend = FALSE) # No legend needed with direct colors
          )
        return(p_var)
      }
    })

    # Observer for cutoff changes to update variability plot colors via proxy
    observe({      
      # Only run for Expression Variability mode and if data cache and inputs are ready
      req(input$y_axis_var == 2, 
          !is.null(variability_plot_data_cache()), 
          !is.null(input$x_cutoff), 
          !is.null(input$y_cutoff))

      cached_data <- variability_plot_data_cache()
      current_x_cutoff <- input$x_cutoff
      current_y_cutoff <- input$y_cutoff

      # Recalculate colors based on new cutoffs
      new_colors <- ifelse(
        cached_data$Variability >= current_y_cutoff & cached_data$log2_Avg_Expression_plus_1 >= current_x_cutoff,
        "blue",
        "grey"
      )
      
      # Update selected_variable_genes reactiveVal
      updated_selected_genes <- dplyr::filter(cached_data, Variability >= current_y_cutoff, log2_Avg_Expression_plus_1 >= current_x_cutoff)
      selected_variable_genes(updated_selected_genes)
      
      log_info(paste("Proxy update: x_cutoff=", current_x_cutoff, ", y_cutoff=", current_y_cutoff, ", num_selected=", nrow(updated_selected_genes)))

      plotlyProxy("variable_genes_plot", session) %>%
        plotlyProxyInvoke("restyle", list(marker = list(color = list(new_colors))))
    })

    output$number_of_genes <- renderText({
      num_genes <- nrow(selected_variable_genes())
      paste0("Number of Selected Genes: ", ifelse(is.null(num_genes), 0, num_genes))
    })

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
      # Groups are not needed for the new variability metric
      # groups <- dataset$selected_groups() 
      # options <- c()
      # for (group in groups) {
      #   options <- c(options, group)
      # }
      
      output$y_axis_operations <- renderUI({
        if(input$y_axis_var == 2) { # This is now the "Expression Variability" option
          # No Group 1 or Group 2 selection needed
          NULL 
        } else {
          # ANOVA still doesn't need group selection here for its specific controls
          NULL
        }
      })
      
      output$y_cutoff <- renderUI({
        if(input$y_axis_var == 1) { # ANOVA
          tagList(
            numericInput(NS(id, "y_cutoff"), "P-value Cutoff", value = 0.05, min = 0, max = 1, step = 0.05),
            actionButton(NS(id, "run_analysis"), "Run Analysis") # Run analysis button only for ANOVA for now
          )
        } else { # Expression Variability
          # For variability, analysis/calculation happens directly in plot rendering based on cutoffs
          # No separate "Run Analysis" button here, plot updates reactively to cutoffs.
          numericInput(NS(id, "y_cutoff"), "Variability Cutoff", value = 1, step = 0.1,min = 0)
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
                 selectInput(NS(id, "y_axis_var"), "Variable Gene Selection Metric",
                             c("ANOVA (p-value)" = 1, "LFC" = 2), selected = 1),
                 uiOutput(NS(id, "y_axis_operations")),
                 numericInput(NS(id, "x_cutoff"), "log2(Avg Expression+1) Cutoff", value = 1, min = 0, step = 0.1),
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
