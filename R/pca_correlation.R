library(ggplot2)
library(corrr)
library(plotly)
library(ggcorrplot)
library(reshape2)
library(shinybusy)
library(heatmaply)
library(logger)

PCACorrelationTabServer <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    log_info("PCA and Correlation tab initialized")

    corr_download_data <- reactiveVal()
    pca_download_data <- reactiveVal()

    get_pca_data <- function(df, groups, selected_samples) {
      log_info("Processing PCA data with {length(selected_samples)} samples")

      df <- df %>%
        select(c('Symbols', 'Genes', all_of(selected_samples)))

      # only use the genes with row sum geater than 0
      rs <- rowSums(df[, 3:ncol(df)])
      use <- (rs > 0)
      df <- df[use,]

      log_info("Filtered to {nrow(df)} genes with non-zero expression")

      if (input$pca_grouped) {
        gene_data_long <- pivot_longer(df, cols = -names(df)[1:2], names_to = "Samples", values_to = "Values")
        gene_data_long_grouped <- gene_data_long %>%
          left_join(groups, by = "Samples")

        gene_data_long_grouped_sum <- gene_data_long_grouped %>%
          group_by(Genes, Symbols, EXPERIMENTAL_GROUP) %>%
          summarize(Values = sum(Values), .groups = "drop")

        result_data <- pivot_wider(gene_data_long_grouped_sum, names_from = EXPERIMENTAL_GROUP, values_from = Values, values_fill = list(value = 0))

        df_transposed <- t(result_data[, 3:ncol(result_data)])

      } else {
        df_transposed <- t(df[, 3:ncol(df)])
      }

      return(df_transposed)
    }

    create_pca_plot <- function() {
      log_info("Creating PCA plot: PC{input$pc_x} vs PC{input$pc_y}")

      req(input$pc_x, input$pc_y, dataset$expression_data(), dataset$groups_data(), dataset$selected_samples_data())
      df <- dataset$filtered_data()
      groups <- dataset$groups_data()
      selected_samples <- dataset$selected_samples_data()

      colnames(df)[1] <- "Symbols"
      colnames(df)[2] <- "Genes"
      colnames(groups)[1] <- "Samples"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"

      pc_x <- as.numeric(input$pc_x)
      pc_y <- as.numeric(input$pc_y)

      df_transposed <- get_pca_data(df, groups, selected_samples)
      pca_comp <- prcomp(df_transposed, scale. = TRUE, center = TRUE)
      percentVar <- pca_comp$sdev^2 / sum(pca_comp$sdev^2)
      pca_df <- data.frame(pca_comp$x, sampleLab = rownames(pca_comp$x), check.names = FALSE)

      if (input$pca_grouped) {
        group_color <- groups %>%
          group_by(Color)
        pca_df <- pca_df %>%
          left_join(group_color, by = c("sampleLab" = "EXPERIMENTAL_GROUP"))
        colors <- unique(pca_df$Color)
        names(colors) <- unique(pca_df$EXPERIMENTAL_GROUP)
        scale_color <- scale_color_manual(values = colors)
      } else {
        pca_df <- pca_df %>%
          left_join(groups, by = c("sampleLab" = "Samples"))
        colors <- pca_df$Color
        names(colors) <- pca_df$sampleLab
        scale_color <- scale_color_manual(values = colors)
        #scale_color <- scale_color_viridis(discrete = TRUE)
      }
      pca_download_data(pca_df)
      p <- pca_df %>%
        ggplot(aes(x = pca_df[, pc_x], y = pca_df[, pc_y], label = sampleLab, color = sampleLab,
                   text = paste("Sample:", sampleLab, "<br>PC", pc_x, ":", pca_df[, pc_x], "<br>PC", pc_y, ":", pca_df[, pc_y]))) +
        geom_point(size = 5) +
        labs(x = paste0("PC", pc_x, ": ", round(percentVar[pc_x] * 100, 1), "% variance"),
             y = paste0("PC", pc_y, ": ", round(percentVar[pc_y] * 100, 1), "% variance")) +
        scale_color +
        theme_minimal() +
        theme(legend.position = "right") +
        theme_linedraw(base_size = 16) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(size = 10),
              legend.position = "right")

      ggplotly(p, tooltip = "text")
      log_info("PCA plot created successfully")
    }

    create_correlation_plot <- function(data, groups, Pearson_or_Spearman) {
      log_info("Creating correlation plot using {Pearson_or_Spearman} coefficient")
      numerical_data <- data[, 3:ncol(data)]
      if (Pearson_or_Spearman == "pearson") {
        title_hierarchical <- "Pearson's Correlation Matrix"
        cor_matrix <- cor(numerical_data, method = "pearson")
      } else {
        title_hierarchical <- "Spearman's Correlation Matrix"
        cor_matrix <- cor(numerical_data, method = "spearman")
      }
      breaks <- seq(0, 1, length.out = 100)

      if (input$pca_grouped) {
        annotation_df <- data.frame(Samples = groups["EXPERIMENTAL_GROUP"], color = groups["Color"])
        annotation_df <- unique(annotation_df);
        annotation_row <- data.frame(Samples = annotation_df$EXPERIMENTAL_GROUP)
        rownames(annotation_row) <- annotation_df$EXPERIMENTAL_GROUP
        annotation_colors <- list(
          Samples = setNames(annotation_df$Color, annotation_df$EXPERIMENTAL_GROUP)
        )
        filtered_row_side_colors <- rownames(cor_matrix)
      } else {
        annotation_df <- data.frame(Samples = groups["Samples"], color = groups["Color"])
        annotation_row <- data.frame(Samples = annotation_df$Samples)
        rownames(annotation_row) <- annotation_df$Samples
        annotation_colors <- list(
          Samples = setNames(annotation_df$Color, annotation_df$Samples)
        )
        filtered_row_side_colors <- annotation_row$Samples[annotation_row$Samples %in% rownames(cor_matrix)]
      }

      corr_download_data(cor_matrix)
      custom_colors <- colorRampPalette(c("#5074AF", "#FFFFFF", "#FFFF66", "#CA4938"))(100)
      p <- heatmaply(
        cor_matrix,
        Rowv = FALSE,
        Colv = FALSE,
        row_side_colors = rownames(cor_matrix),
        row_side_palette = annotation_colors$Samples,
        main = title_hierarchical,
        colors = custom_colors(100),
        scale_fill_gradient_fun = ggplot2::scale_fill_gradientn(colors = custom_colors, limits = c(input$scale[1], input$scale[2])),
        fontsize_row = 10,
        fontsize_col = 10,
        grid_color = "grey",
        xlab = NULL,
        ylab = NULL,
        colorbar_thickness = 10,
        label_names = c("X-axis", "Y-axis", "Value")
      ) %>% layout(showlegend = FALSE)

      log_info("Correlation plot created successfully")
      return(p)
    }


    create_clustering_plot <- function(data, groups, Pearson_or_Spearman) {
      data_sample_expr <- data[, -c(1, 2)]

      if (Pearson_or_Spearman == "pearson") {
        title_kmeans <- "K-means Clustered, Pearson's Correlation Matrix"
        doc_label_kmeans <- "Kmeans_Pearsons_Correlation_Matrix"
        title_hierarchical <- "Hierarchically Clustered, Pearson's Correlation Matrix"
        doc_label_hierarchical <- "Hierarchical_Pearsons_Correlation_Matrix"
        cluster_assignment <- "Kmeans_Pearsons_Correlation_Matrix_cluster_assignments"
        cor_matrix <- cor(data_sample_expr, method = "pearson")
      } else {
        title_kmeans <- "K-means Clustered, Spearman's Correlation Matrix"
        doc_label_kmeans <- "Kmeans_Spearmans_Correlation_Matrix"
        title_hierarchical <- "Hierarchically Clustered, Spearman's Correlation Matrix"
        doc_label_hierarchical <- "Hierarchical_Spearmans_Correlation_Matrix"
        cluster_assignment <- "Kmeans_Spearmans_Correlation_Matrix_cluster_assignments"
        cor_matrix <- cor(data_sample_expr, method = "spearman")
      }

      set.seed(40)
      breaks <- seq(0, 1, length.out = 100)

      if (input$pca_grouped) {
        annotation_df <- data.frame(Samples = groups["EXPERIMENTAL_GROUP"], color = groups["Color"])
        annotation_df <- unique(annotation_df);
        annotation_row <- data.frame(Samples = annotation_df$EXPERIMENTAL_GROUP)
        rownames(annotation_row) <- annotation_df$EXPERIMENTAL_GROUP
        annotation_colors <- list(
          Samples = setNames(annotation_df$Color, annotation_df$EXPERIMENTAL_GROUP)
        )
      } else {
        annotation_df <- data.frame(Samples = groups["Samples"], color = groups["Color"])
        annotation_row <- data.frame(Samples = annotation_df$Samples)
        rownames(annotation_row) <- annotation_df$Samples
        annotation_colors <- list(
          Samples = setNames(annotation_df$Color, annotation_df$Samples)
        )
      }

      distance_matrix <- function(cor_matrix) {
        as.dist(1 - cor_matrix)
      }

      corr_download_data(distance_matrix(cor_matrix))
      custom_colors <- colorRampPalette(c("#5074AF", "#FFFFFF", "#FFFF66", "#CA4938"))(100)
      p <- heatmaply(
        cor_matrix,
        row_side_colors = rownames(cor_matrix),
        row_side_palette = annotation_colors$Samples,
        Rowv = hclust(distance_matrix(cor_matrix)),
        Colv = hclust(distance_matrix(t(cor_matrix))),
        main = title_hierarchical,
        colors = custom_colors(100),
        scale_fill_gradient_fun = ggplot2::scale_fill_gradientn(colors = custom_colors, limits = c(input$scale[1], input$scale[2])),
        fontsize_row = 10,
        fontsize_col = 10,
        grid_color = "grey",
        xlab = NULL,
        ylab = NULL,
        colorbar_thickness = 10,
        label_names = c("X-axis", "Y-axis", "Value")
      ) %>% layout(showlegend = FALSE)
      return(p)
    }


    output$pca <- renderPlotly({
      req(dataset$expression_data())
      create_pca_plot()
    })

    output$correlation <- renderPlotly({
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

      if (input$clustered) {
        plot <- create_clustering_plot(data, groups, input$coefficient)
      } else {
        plot <- create_correlation_plot(data, groups, input$coefficient)
      }

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
        width: 100%; /* 容器宽度占满父元素宽度 */
        height: 0; /* 初始高度设为0 */
        padding-bottom: 100%; /* 高度等于宽度的百分比 */
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
      }
    "))
    ),
    titlePanel("PCA and Correlation"),
    fluidRow(
      column(6,
             fluidRow(class = "bottom-centered",
                      column(2,
                             selectInput(NS(id, "pc_x"), "X Axis:",
                                         c("PC1" = 1), selected = 1)
                      ),
                      column(2,
                             selectInput(NS(id, "pc_y"), "Y Axis:",
                                         c("PC2" = 2), selected = 2)
                      ),
                      column(2,
                             materialSwitch(inputId = NS(id, "pca_grouped"), label = "Show by group: ", value = FALSE, status = "primary")
                      ),
                      column(3,
                             downloadButton(NS(id, "download_pca"), "Download Data")
                      ),
             ),
             plotlyOutput(NS(id, "pca"), width = "100%", height = "700px")
      ),
      column(6,
             fluidRow(class = "bottom-centered",
                      column(4,
                             selectInput(NS(id, "coefficient"), "Coefficient:",
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
                      column(3,
                             materialSwitch(inputId = NS(id, "clustered"), label = "Hierarchical Clustering: ", value = FALSE, status = "primary")
                      ),
                      column(2,
                             downloadButton(NS(id, "download_corr"), "Download Data")
                      )
             ),
             add_busy_spinner(spin = "fading-circle", color = "#000000"),
             div(id = "plot-container",
                 plotlyOutput(NS(id, "correlation"), height = "85%")
             )

      )
    )
  )
}

