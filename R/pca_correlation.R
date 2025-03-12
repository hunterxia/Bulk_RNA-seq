library(ggplot2)
library(corrr)
library(plotly)
library(ggcorrplot)
library(reshape2)
library(shinybusy)
library(heatmaply)
library(logger)
library(pheatmap)

# Define server-side module for the PCA and Correlation tab
PCACorrelationTabServer <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    log_info("PCA and Correlation tab initialized")

    # Reactive values to store data for downloads and heatmap object
    corr_download_data <- reactiveVal()
    pca_download_data <- reactiveVal()
    heatmap_obj <- reactiveVal()

    # Function to get PCA data
    get_pca_data <- function(df, groups, selected_samples) {
      log_info("Processing PCA data with {length(selected_samples)} samples")

      # select only the selected samples
      df <- df %>%
        select(c('Symbols', 'Genes', all_of(selected_samples)))

      # only use the genes with row sum geater than 0
      rs <- rowSums(df[, 3:ncol(df)])
      use <- (rs > 0)
      df <- df[use,]

      log_info("Filtered to {nrow(df)} genes with non-zero expression")

      # If PCA grouping is enabled, transform the data by aggregating by experimental group
      if (input$pca_grouped) {
        gene_data_long <- pivot_longer(df, cols = -names(df)[1:2], names_to = "Samples", values_to = "Values")
        # Merge with group information
        gene_data_long_grouped <- gene_data_long %>%
          left_join(groups, by = "Samples")

        # Merge with group information
        gene_data_long_grouped_sum <- gene_data_long_grouped %>%
          group_by(Genes, Symbols, EXPERIMENTAL_GROUP) %>%
          summarize(Values = sum(Values), .groups = "drop")

        # Reshape data back to wide format with groups as columns
        result_data <- pivot_wider(gene_data_long_grouped_sum, names_from = EXPERIMENTAL_GROUP, values_from = Values, values_fill = list(value = 0))

        # Transpose the data for PCA (samples as rows, genes as columns)
        df_transposed <- t(result_data[, 3:ncol(result_data)])

      } else {
        df_transposed <- t(df[, 3:ncol(df)])
      }

      return(df_transposed)
    }

    # create_pca_plot <- function() {
    #   log_info("Creating PCA plot: PC{input$pc_x} vs PC{input$pc_y}")
    #
    #   req(input$pc_x, input$pc_y, dataset$expression_data(), dataset$groups_data(), dataset$selected_samples_data())
    #   df <- dataset$filtered_data()
    #   groups <- dataset$groups_data()
    #   selected_samples <- dataset$selected_samples_data()
    #
    #   colnames(df)[1] <- "Symbols"
    #   colnames(df)[2] <- "Genes"
    #   colnames(groups)[1] <- "Samples"
    #   colnames(groups)[2] <- "EXPERIMENTAL_GROUP"
    #
    #   pc_x <- as.numeric(input$pc_x)
    #   pc_y <- as.numeric(input$pc_y)
    #
    #   df_transposed <- get_pca_data(df, groups, selected_samples)
    #   pca_comp <- prcomp(df_transposed, scale. = TRUE, center = TRUE)
    #   percentVar <- pca_comp$sdev^2 / sum(pca_comp$sdev^2)
    #   pca_df <- data.frame(pca_comp$x, sampleLab = rownames(pca_comp$x), check.names = FALSE)
    #
    #   if (input$pca_grouped) {
    #     group_color <- groups %>%
    #       group_by(Color)
    #     pca_df <- pca_df %>%
    #       left_join(group_color, by = c("sampleLab" = "EXPERIMENTAL_GROUP"))
    #     colors <- unique(pca_df$Color)
    #     names(colors) <- unique(pca_df$EXPERIMENTAL_GROUP)
    #     scale_color <- scale_color_manual(values = colors)
    #   } else {
    #     pca_df <- pca_df %>%
    #       left_join(groups, by = c("sampleLab" = "Samples"))
    #     colors <- pca_df$Color
    #     names(colors) <- pca_df$sampleLab
    #     scale_color <- scale_color_manual(values = colors)
    #     #scale_color <- scale_color_viridis(discrete = TRUE)
    #   }
    #   pca_download_data(pca_df)
    #
    #   x_col <- paste0("PC", pc_x)
    #   y_col <- paste0("PC", pc_y)
    #
    #   plotly_obj <- plot_ly(
    #     data = pca_df,
    #     x = ~get(x_col),
    #     y = ~get(y_col),
    #     type = 'scatter',
    #     mode = 'markers',
    #     color = ~color_var,
    #     colors = color_mapping,
    #     text = ~paste("Sample:", sampleLab,
    #                   "<br>", x_col, ":", round(get(x_col), 2),
    #                   "<br>", y_col, ":", round(get(y_col), 2)),
    #     hoverinfo = "text"
    #   ) %>%
    #     layout(
    #       title = sprintf("PCA Plot: %s vs %s", x_col, y_col),
    #       xaxis = list(title = paste0(x_col, ": ", round(percentVar[pc_x] * 100, 1), "% variance")),
    #       yaxis = list(title = paste0(y_col, ": ", round(percentVar[pc_y] * 100, 1), "% variance")),
    #       legend = list(orientation = "v", x = 1, y = 1)
    #     )
    #
    #   log_info("PCA plot created successfully")
    #
    #   plotly_obj


    # p <- pca_df %>%
    #   ggplot(aes(x = pca_df[, pc_x], y = pca_df[, pc_y], label = sampleLab, color = sampleLab,
    #              text = paste("Sample:", sampleLab, "<br>PC", pc_x, ":", pca_df[, pc_x], "<br>PC", pc_y, ":", pca_df[, pc_y]))) +
    #   geom_point(size = 5) +
    #   labs(x = paste0("PC", pc_x, ": ", round(percentVar[pc_x] * 100, 1), "% variance"),
    #        y = paste0("PC", pc_y, ": ", round(percentVar[pc_y] * 100, 1), "% variance")) +
    #   scale_color +
    #   theme_minimal() +
    #   theme(legend.position = "right") +
    #   theme_linedraw(base_size = 16) +
    #   theme(panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank(),
    #         legend.title = element_blank(),
    #         legend.text = element_text(size = 10),
    #         legend.position = "right")
    #
    # ggplotly(p, tooltip = "text")
    # log_info("PCA plot created successfully")
    # }

    # Function to create the PCA plot using plotly
    create_pca_plot <- function() {
      log_info(sprintf("Creating PCA plot: PC%s vs PC%s", input$pc_x, input$pc_y))

      # Require that necessary inputs and data are available
      req(input$pc_x, input$pc_y,
          dataset$expression_data(),
          dataset$groups_data(),
          dataset$selected_samples_data())

      # Retrieve filtered expression data and group info
      df <- dataset$filtered_data()
      groups <- dataset$groups_data()
      selected_samples <- dataset$selected_samples_data()

      # Rename columns for consistency
      colnames(df)[1:2] <- c("Symbols", "Genes")
      colnames(groups)[1:2] <- c("Samples", "EXPERIMENTAL_GROUP")

      # Convert selected PC indices to numeric
      pc_x <- as.numeric(input$pc_x)
      pc_y <- as.numeric(input$pc_y)

      # Get transposed data for PCA calculation
      df_transposed <- get_pca_data(df, groups, selected_samples)

      # Perform Principal Component Analysis with scaling and centering
      pca_comp <- prcomp(df_transposed, scale. = TRUE, center = TRUE)
      percentVar <- pca_comp$sdev^2 / sum(pca_comp$sdev^2)

      # Create a data frame with PCA components and sample labels
      pca_df <- data.frame(pca_comp$x, sampleLab = rownames(pca_comp$x), check.names = FALSE)

      # Depending on whether grouping is enabled, merge group colors and set up color mapping
      if (input$pca_grouped) {
        pca_df <- pca_df %>% left_join(groups, by = c("sampleLab" = "EXPERIMENTAL_GROUP"))

        if ("Color" %in% colnames(groups)) {
          group_color_map <- groups %>%
            select(EXPERIMENTAL_GROUP, Color) %>%
            distinct()
          colors <- group_color_map$Color
          names(colors) <- group_color_map$EXPERIMENTAL_GROUP
          color_mapping <- colors
        } else {
          color_mapping <- "Viridis"
        }

        color_var <- pca_df$sampleLab
      } else {
        pca_df <- pca_df %>% left_join(groups, by = c("sampleLab" = "Samples"))
        if ("Color" %in% colnames(pca_df)) {
          sample_color_map <- pca_df %>%
            select(sampleLab, Color) %>%
            distinct()
          colors <- sample_color_map$Color
          names(colors) <- sample_color_map$sampleLab
          if (length(colors) < length(unique(pca_df$sampleLab))) {
            unique_samples <- unique(pca_df$sampleLab)
            colors <- viridis::viridis(length(unique_samples))
            names(colors) <- unique_samples
          }
          color_mapping <- colors
        } else {
          color_mapping <- "Viridis"
        }
        color_var <- pca_df$sampleLab
      }

      # Save PCA data for download
      pca_download_data(pca_df)

      # Determine which principal components to plot
      x_col <- paste0("PC", pc_x)
      y_col <- paste0("PC", pc_y)

      # Create a plotly scatter plot for PCA
      plotly_obj <- plot_ly(
        data = pca_df,
        x = ~get(x_col),
        y = ~get(y_col),
        type = 'scatter',
        mode = 'markers',
        marker = list(size = 15),
        color = ~color_var,
        colors = color_mapping,
        text = ~paste("Sample:", sampleLab,
                      "<br>", x_col, ":", round(get(x_col), 2),
                      "<br>", y_col, ":", round(get(y_col), 2)),
        hoverinfo = "text"
      ) %>%
        layout(
          title = sprintf("PCA Plot: %s vs %s", x_col, y_col),
          xaxis = list(title = paste0(x_col, ": ", round(percentVar[pc_x] * 100, 1), "% variance")),
          yaxis = list(title = paste0(y_col, ": ", round(percentVar[pc_y] * 100, 1), "% variance")),
          legend = list(orientation = "v", x = 1, y = 1)
        )

      log_info("PCA plot created successfully")

      plotly_obj
    }

    # create_correlation_plot <- function(data, groups, Pearson_or_Spearman) {
    #   log_info("Creating correlation plot using {Pearson_or_Spearman} coefficient")
    #   numerical_data <- data[, 3:ncol(data)]
    #
    #   if (Pearson_or_Spearman == "pearson") {
    #     title_hierarchical <- "Pearson's Correlation Matrix"
    #     cor_matrix <- cor(numerical_data, method = "pearson")
    #   } else {
    #     title_hierarchical <- "Spearman's Correlation Matrix"
    #     cor_matrix <- cor(numerical_data, method = "spearman")
    #   }
    #   breaks <- seq(0, 1, length.out = 100)
    #
    #   if (input$pca_grouped) {
    #     annotation_df <- data.frame(Samples = groups["EXPERIMENTAL_GROUP"], color = groups["Color"])
    #     annotation_df <- unique(annotation_df);
    #     annotation_row <- data.frame(Samples = annotation_df$EXPERIMENTAL_GROUP)
    #     rownames(annotation_row) <- annotation_df$EXPERIMENTAL_GROUP
    #     annotation_colors <- list(
    #       Samples = setNames(annotation_df$Color, annotation_df$EXPERIMENTAL_GROUP)
    #     )
    #     filtered_row_side_colors <- rownames(cor_matrix)
    #   } else {
    #     annotation_df <- data.frame(Samples = groups["Samples"], color = groups["Color"])
    #     annotation_row <- data.frame(Samples = annotation_df$Samples)
    #     rownames(annotation_row) <- annotation_df$Samples
    #     annotation_colors <- list(
    #       Samples = setNames(annotation_df$Color, annotation_df$Samples)
    #     )
    #     filtered_row_side_colors <- annotation_row$Samples[annotation_row$Samples %in% rownames(cor_matrix)]
    #   }
    #
    #   corr_download_data(cor_matrix)
    #   custom_colors <- colorRampPalette(c("#5074AF", "#FFFFFF", "#FFFF66", "#CA4938"))(100)
    #   p <- heatmaply(
    #     cor_matrix,
    #     Rowv = FALSE,
    #     Colv = FALSE,
    #     row_side_colors = rownames(cor_matrix),
    #     row_side_palette = annotation_colors$Samples,
    #     main = title_hierarchical,
    #     colors = custom_colors(100),
    #     scale_fill_gradient_fun = ggplot2::scale_fill_gradientn(colors = custom_colors, limits = c(input$scale[1], input$scale[2])),
    #     fontsize_row = 10,
    #     fontsize_col = 10,
    #     grid_color = "grey",
    #     xlab = NULL,
    #     ylab = NULL,
    #     colorbar_thickness = 10,
    #     label_names = c("X-axis", "Y-axis", "Value")
    #   ) %>% layout(showlegend = FALSE)
    #
    #   log_info("Correlation plot created successfully")
    #   return(p)
    # }

    # Function to create the correlation heatmap using pheatmap
    create_correlation_plot <- function(data, groups, Pearson_or_Spearman) {
      log_info("Creating correlation plot using {Pearson_or_Spearman} coefficient")

      numerical_data <- data[, 3:ncol(data)]

      # Compute the correlation matrix based on the specified method
      if (Pearson_or_Spearman == "pearson") {
        title_hierarchical <- "Pearson's Correlation Matrix"
        cor_matrix <- cor(numerical_data, method = "pearson", use = "pairwise.complete.obs")
      } else {
        title_hierarchical <- "Spearman's Correlation Matrix"
        cor_matrix <- cor(numerical_data, method = "spearman", use = "pairwise.complete.obs")
      }

      # Create a custom color palette for the heatmap
      custom_color_fun <- colorRampPalette(c("#5074AF", "#FFFFFF", "#FFFF66", "#CA4938"))
      custom_colors <- custom_color_fun(100)

      # Create a custom color palette for the heatmap
      scale_limits <- if (!is.null(input$scale)) {
        c(input$scale[1], input$scale[2])
      } else {
        c(min(cor_matrix, na.rm = TRUE), max(cor_matrix, na.rm = TRUE))
      }
      breaks <- seq(scale_limits[1], scale_limits[2], length.out = 101)

      # Set up annotation data and colors based on grouping option
      if (input$pca_grouped) {
        # Create annotation data frame from group information for experimental groups
        annotation_df <- data.frame(
          Samples = groups$EXPERIMENTAL_GROUP,
          Color = groups$Color,
          stringsAsFactors = FALSE
        )
        annotation_df <- unique(annotation_df)

        annotation_row <- data.frame(Samples = annotation_df$Samples, row.names = annotation_df$Samples, stringsAsFactors = FALSE)

        row_annotation <- annotation_row[rownames(annotation_row) %in% rownames(cor_matrix), , drop = FALSE]

        filtered_samples <- rownames(row_annotation)

        annotation_colors <- list(
          Samples = setNames(annotation_df$Color, annotation_df$Samples)[filtered_samples]
        )
      } else {
        annotation_df <- data.frame(
          Samples = groups$Samples,
          Color = groups$Color,
          stringsAsFactors = FALSE
        )
        annotation_df <- unique(annotation_df)
        annotation_row <- data.frame(Samples = annotation_df$Samples, row.names = annotation_df$Samples, stringsAsFactors = FALSE)
        row_annotation <- annotation_row[rownames(annotation_row) %in% rownames(cor_matrix), , drop = FALSE]
        filtered_samples <- rownames(row_annotation)
        annotation_colors <- list(
          Samples = setNames(annotation_df$Color, annotation_df$Samples)[filtered_samples]
        )
      }

      # Save the correlation matrix for download purposes
      corr_download_data(cor_matrix)

      # Create the heatmap using pheatmap with silent = TRUE to capture the object
      pheatmap_obj <- pheatmap::pheatmap(
        mat = cor_matrix,
        color = custom_colors,
        breaks = breaks,
        main = title_hierarchical,
        annotation_row = row_annotation,
        annotation_colors = annotation_colors,
        cluster_rows = input$hierarchical_clustering,
        cluster_cols = input$hierarchical_clustering,
        fontsize_row = 10,
        fontsize_col = 10,
        border_color = "grey",
        silent = TRUE
      )

      # Store the heatmap object for download
      heatmap_obj(pheatmap_obj)
      log_info("Correlation plot created successfully")

      return(pheatmap_obj)
    }

    # Render the PCA plot using Plotly
    output$pca <- renderPlotly({
      req(dataset$expression_data())
      create_pca_plot()
    })

    # Render the correlation heatmap plot
    output$correlation <- renderPlot({
      req(dataset$filtered_data(), dataset$groups_data())
      gene_data <- dataset$filtered_data()
      groups <- dataset$groups_data()
      colnames(gene_data)[1] <- "Genes"
      colnames(gene_data)[2] <- "Symbols"
      colnames(groups)[1] <- "Samples"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"

      # show by grouped
      if (input$pca_grouped) {
        gene_data_long <- pivot_longer(gene_data, cols = -names(gene_data)[1:2], names_to = "Samples", values_to = "Values")

        gene_data_long_grouped <- gene_data_long %>%
          left_join(groups, by = "Samples")

        gene_data_long_grouped_sum <- gene_data_long_grouped %>%
          group_by(Genes, Symbols, EXPERIMENTAL_GROUP) %>%
          summarize(Values = sum(Values), .groups = "drop")

        data <- pivot_wider(gene_data_long_grouped_sum, names_from = EXPERIMENTAL_GROUP, values_from = Values, values_fill = list(value = 0))
      } else {
        data <- gene_data
      }

      # if (input$clustered) {
      #   plot <- create_clustering_plot(data, groups, input$coefficient)
      # } else {
      #   plot <- create_correlation_plot(data, groups, input$coefficient)
      # }
      plot <- create_correlation_plot(data, groups, input$coefficient)

      return(plot)

    })

    # download pca data
    output$download_pca <- downloadHandler(
      filename = function() {
        log_info("Downloading PCA data")
        paste("pca-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        req(input$pc_x, input$pc_y)
        write.csv(pca_download_data(), file)
        log_info("PCA data downloaded successfully")
      }
    )

    # download correlation data
    output$download_corr <- downloadHandler(
      filename = function() {
        paste("correlation-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(corr_download_data(), file)
      }
    )

    # download heatmap
    output$download_heatmap <- downloadHandler(
      filename = function() {
        paste("heatmap-", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = 10, height = 8)
        grid::grid.newpage()
        grid::grid.draw(heatmap_obj()$gtable)
        dev.off()
      }
    )

    #update selector options
    observeEvent(dataset$selected_samples(), {
      tryCatch({
        df <- dataset$expression_data()
        groups <- dataset$groups_data()
        selected_samples <- dataset$selected_samples_data()
        colnames(df)[1] <- "Genes"
        colnames(df)[2] <- "Symbols"
        colnames(groups)[1] <- "Samples"
        colnames(groups)[2] <- "EXPERIMENTAL_GROUP"
        df_transposed <- get_pca_data(df, groups, selected_samples)
        pca_comp <- prcomp(df_transposed, scale. = TRUE, center = TRUE)
        col_names <- colnames(pca_comp$x)
        options <- c()
        for (idx in seq_along(col_names)) {
          col_name <- col_names[idx]
          options[col_name] = idx
        }
        updateSelectInput(session, "pc_x", choices = options, selected = 1)
        updateSelectInput(session, "pc_y", choices = options, selected = 2)
        log_info("PCA selector options updated with {length(options)} choices")
      }, error = function(e) {
        log_error("Error updating PCA selector options: {e$message}")
      })

    })

    #update selector options
    observeEvent(input$pca_grouped, {
      req(dataset$expression_data(), dataset$groups_data(), dataset$selected_samples_data())
      df <- dataset$expression_data()
      groups <- dataset$groups_data()
      selected_samples <- dataset$selected_samples_data()
      colnames(df)[1] <- "Genes"
      colnames(df)[2] <- "Symbols"
      colnames(groups)[1] <- "Samples"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"
      df_transposed <- get_pca_data(df, groups, selected_samples)
      pca_comp <- prcomp(df_transposed, scale. = TRUE, center = TRUE)
      col_names <- colnames(pca_comp$x)
      options <- c()
      for (idx in seq_along(col_names)) {
        col_name <- col_names[idx]
        options[col_name] = idx
      }
      updateSelectInput(session, "pc_x", choices = options, selected = 1)
      updateSelectInput(session, "pc_y", choices = options, selected = 2)
    })

  })
}

PCACorrelationTabUI <- function(id) {
  fluidPage(
    tags$head(
      tags$style(HTML("
      #plot-container {
        width: 100%;
        height: 0; /* 初始高度设为0 */
        padding-bottom: 100%;
        position: relative;
      }
      #pca_correlation-correlation {
        position: absolute;
        top: 0;
        bottom: 0;
        left: 0;
        right: 0;
      }
      .bottom-centered {
        display: flex;
        align-items: center;
        flex-direction: row;
      }
      .inline-row {
        display: flex;
        flex-direction: row;
        align-items: center;
        justify-content: center;

      }


    "))
    ),
    titlePanel("PCA and Correlation"),
    # fluidRow(
    #   column(6,
    #          fluidRow(class = "bottom-centered",
    #                   column(3,
    #                          materialSwitch(inputId = NS(id, "pca_grouped"), label = "Show by group: ", value = FALSE, status = "primary")
    #                   ),
    #                   column(3,
    #                          selectInput(NS(id, "pc_x"), "X Axis:",
    #                                      c("PC1" = 1), selected = 1)
    #                   ),
    #                   column(3,
    #                          selectInput(NS(id, "pc_y"), "Y Axis:",
    #                                      c("PC2" = 2), selected = 2)
    #                   ),
    #                   column(3,
    #                          downloadButton(NS(id, "download_pca"), "Download PC coordinates")
    #                   ),
    #          ),
    #          plotlyOutput(NS(id, "pca"), width = "100%", height = "700px")
    #   ),
    #   column(6,
    #          fluidRow(class = "bottom-centered",
    #                   column(4,
    #                          selectInput(NS(id, "coefficient"), "Coefficient:",
    #                                      c("Pearson’s coefficient" = "pearson",
    #                                        " Spearman’s rank coefficient" = "spearman"))
    #                   ),
    #                   column(3,
    #                          sliderInput(NS(id, "scale"),
    #                                      "Scale:",
    #                                      min = -1,
    #                                      max = 1,
    #                                      step = 0.05,
    #                                      value = c(0, 1))
    #                   ),
    #                   # column(3,
    #                   #        materialSwitch(inputId = NS(id, "clustered"), label = "Hierarchical Clustering: ", value = FALSE, status = "primary")
    #                   # ),
    #                   column(2,
    #                          downloadButton(NS(id, "download_corr"), "Download correlations")
    #                   )
    #          ),
    #          add_busy_spinner(spin = "fading-circle", color = "#000000"),
    #          div(id = "plot-container",
    #              plotOutput(NS(id, "correlation"), height = "85%")
    #          )
    #
    #   )
    # )

    fluidRow(
      column(6,
             fluidRow(
               column(6,
                      materialSwitch(inputId = NS(id, "pca_grouped"),
                                     label = "Show by group: ",
                                     value = FALSE,
                                     status = "primary")
               ),
               column(6, materialSwitch(inputId = NS(id, "hierarchical_clustering"),
                                        label = "Hierarchical Clustering: ",
                                        value = FALSE,
                                        status = "primary")
               )
             ),
             fluidRow(
               column(12,
                      downloadButton(NS(id, "download_heatmap"), "Download Heatmap")
               )
             ),
             fluidRow(
               class = "inline-row",
               column(4,
                      selectInput(NS(id, "pc_x"),
                                  "X Axis:",
                                  c("PC1" = 1),
                                  selected = 1)
               ),
               column(4,
                      selectInput(NS(id, "pc_y"),
                                  "Y Axis:",
                                  c("PC2" = 2),
                                  selected = 2)
               ),
               column(4,
                      downloadButton(NS(id, "download_pca"), "Download PC coordinates")
               ),
             ),
             plotlyOutput(NS(id, "pca"), width = "100%", height = "700px")
      ),
      column(6,
             fluidRow(
               class = "inline-row",
               column(4,
                      selectInput(NS(id, "coefficient"),
                                  "Coefficient:",
                                  c("Pearson’s coefficient" = "pearson",
                                    " Spearman’s rank coefficient" = "spearman"))
               ),
               column(3,
                      sliderInput(NS(id, "scale"),
                                  "Scale:",
                                  min = -1,
                                  max = 1,
                                  step = 0.05,
                                  value = c(0, 1))
               ),
               column(2,
                      downloadButton(NS(id, "download_corr"), "Download correlations")
               )
             ),
             add_busy_spinner(spin = "fading-circle", color = "#000000"),
             div(id = "plot-container",
                 plotOutput(NS(id, "correlation"), height = "85%")
             )
      )
    ))
}

