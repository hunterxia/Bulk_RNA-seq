library(ggplot2)
library(corrr)
library(plotly)
library(ggcorrplot)
library(reshape2)
library(shinybusy)

PCACorrelationTabServer <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    get_pca_data <- function(df, groups, selected_samples) {
      df <- df %>%
        select(c('Symbols', 'Genes', all_of(selected_samples)))

      # only use the genes with row sum geater than 0
      rs <- rowSums(df[, 3:ncol(df)])
      use <- (rs > 0)
      df <- df[use, ]
      
      if (input$pca_grouped) {
        gene_data_long <- pivot_longer(df, cols = -names(df)[1:2], names_to = "Samples", values_to = "Values")
        gene_data_long_grouped <- gene_data_long %>%
          left_join(groups, by = "Samples") 
        
        gene_data_long_grouped_sum <- gene_data_long_grouped %>%
          group_by(Genes, Symbols, EXPERIMENTAL_GROUP) %>%
          summarize(Values = sum(Values), .groups = "drop")
        
        result_data <- pivot_wider(gene_data_long_grouped_sum, names_from = EXPERIMENTAL_GROUP, values_from = Values, values_fill = list(value = 0))
        
        df_transposed <- t(result_data[,3:ncol(result_data)])
        
      } else {
        df_transposed <- t(df[,3:ncol(df)])
      }
      
      return (df_transposed)
    }
    
    create_pca_plot <- function(){
      req(input$pc_x, input$pc_y, dataset$expression_data(), dataset$groups_data(), dataset$selected_samples_data())
      df <- dataset$expression_data()
      groups <- dataset$groups_data()
      selected_samples <- dataset$selected_samples_data()
      
      colnames(df)[1] <- "Symbols"
      colnames(df)[2] <- "Genes"
      colnames(groups)[1] <- "Samples"
      colnames(groups)[2] <- "EXPERIMENTAL_GROUP"
      
      pc_x <- as.numeric(input$pc_x)
      pc_y <- as.numeric(input$pc_y)
      # browser()
      df_transposed <- get_pca_data(df, groups, selected_samples)
      pca_comp <- prcomp(df_transposed, scale. = TRUE, center = TRUE)
      percentVar <- pca_comp$sdev^2/sum(pca_comp$sdev^2)
      pca_df <- data.frame(pca_comp$x, sampleLab = rownames(pca_comp$x), check.names = FALSE)
      
      if (input$pca_grouped) {
        group_color <- groups %>%
          group_by(Color) %>%
          slice(1) 
        pca_df <- pca_df %>%
          left_join(group_color, by = c("sampleLab" = "EXPERIMENTAL_GROUP")) 
        colors <- unique(pca_df$Color)
        names(colors) <- unique(pca_df$Var2)
        scale_color <- scale_color_manual(values = colors)
      } else {
        scale_color <- scale_color_viridis(discrete = TRUE)
      }
   
      p <- pca_df %>%
        ggplot(aes(x = pca_df[, pc_x],y = pca_df[, pc_y], label = sampleLab, color = sampleLab, 
                   text = paste("Sample:", sampleLab, "<br>PC", pc_x ,":", pca_df[, pc_x], "<br>PC", pc_y  ,":", pca_df[, pc_y]))) +
        geom_point(size=5) +
        labs(x=paste0("PC", pc_x, ": ", round(percentVar[pc_x]*100,1), "% variance"),
             y=paste0("PC", pc_y, ": ", round(percentVar[pc_y]*100,1), "% variance")) +
        scale_color + 
        theme_minimal() + 
        theme(legend.position = "right") +
        theme_linedraw(base_size=16) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              legend.title = element_blank(),
              legend.text = element_text(size = 10),
              legend.position = "right")
      
      ggplotly(p, tooltip = "text")
    }
    
    create_correlation_plot <- function(data, Pearson_or_Spearman){
      numerical_data <- data[,3:ncol(data)]
      if(Pearson_or_Spearman == "pearson") {
        title_hierarchical <- "Pearson's Correlation Matrix"
        cor_matrix <- cor(numerical_data, method = "pearson")
      } else {
        title_hierarchical <- "Spearman's Correlation Matrix"
        cor_matrix <- cor(numerical_data, method = "spearman")
      }
      p <- pheatmap::pheatmap(cor_matrix, 
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "complete",
               fontsize_row = 10,
               main = title_hierarchical,
               treeheight_row = 30,
               angle_col = 315,
               fontsize_col = 10,
               width = 10,
               height = 10)
      print(p)
      return(p)
    }
    
    
    create_clustering_plot <- function(data, Pearson_or_Spearman) {
      data_sample_expr <- data[,-c(1,2)]

      if(Pearson_or_Spearman == "pearson") {
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
      hierarchical_clustered_correction_plot <- pheatmap::pheatmap(cor_matrix, 
                                                                   cluster_rows = TRUE,
                                                                   cluster_cols = TRUE, 
                                                                   clustering_distance_rows = "correlation",
                                                                   clustering_distance_cols = "correlation",
                                                                   fontsize_row = 10,
                                                                   main = title_hierarchical,
                                                                   treeheight_row = 30,
                                                                   angle_col = 315,
                                                                   fontsize_col = 10,
                                                                   width = 10,
                                                                   height = 10)
      print(hierarchical_clustered_correction_plot)
      return(hierarchical_clustered_correction_plot)
    }
    
    
    output$pca <- renderPlotly({
      req(dataset$expression_data())
      create_pca_plot()
    })
    
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
      
      if (input$clustered) {
        plot <- create_clustering_plot(data, input$coefficient)
      } else {
        plot <- create_correlation_plot(data, input$coefficient)
      }
      
      print(plot)
      return(plot)
      
    }, res = 96)
    
    # download pca data
    output$download_pca <- downloadHandler(
      filename = function() {
        paste("pca-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        req(input$pc_x, input$pc_y)
        pc_x <- as.numeric(input$pc_x)
        pc_y <- as.numeric(input$pc_y)
        df <- dataset$expression_data()
        colnames <- colnames(df)
        df_transposed <- t(df[,3:ncol(df)])
        pca_comp <- prcomp(df_transposed, scale. = TRUE, center = TRUE)
        percentVar <- pca_comp$sdev^2/sum(pca_comp$sdev^2)
        pca_df <- data.frame(pca_comp$x, sampleLab = rownames(pca_comp$x), check.names = FALSE)
        write.csv(pca_df, file)
      }
    )
    
    # download correlation data
    output$download_corr <- downloadHandler(
      filename = function() {
        paste("correlation-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        df <- dataset$expression_data()
        numerical_data <- df[,3:ncol(df)]
        corr_matrix <- cor(data_normalized, method = input$coefficient)
        corr_melted <- melt(corr_matrix)
        write.csv(corr_melted, file)
      }
    )
    
    #update selector options
    observeEvent(dataset$selected_samples(), {
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
        print(options)
      }
      updateSelectInput(session, "pc_x", choices = options, selected=1)
      updateSelectInput(session, "pc_y", choices = options, selected=2)
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
        print(options)
      }
      updateSelectInput(session, "pc_x", choices = options, selected=1)
      updateSelectInput(session, "pc_y", choices = options, selected=2)
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
    "))
    ),
    titlePanel("PCA and Correlation"),
    fluidRow(
      column(6, 
       fluidRow(
         column(2, 
                selectInput(NS(id, "pc_x"), "X Axis:",
                            c("PC1" = 1), selected=1)
         ),
         column(2, 
                selectInput(NS(id, "pc_y"), "Y Axis:",
                            c("PC2" = 2), selected=2)
         ),
         column(2, 
                materialSwitch(inputId = NS(id, "pca_grouped"), label = "Show by group: ", value = FALSE, status = "primary")
         )
       ),
       downloadButton(NS(id, "download_pca"), "Download Data"),
       plotlyOutput(NS(id, "pca"), width = "100%", height = "700px")
      ),
      column(6, 
       fluidRow(
         column(4, 
                selectInput(NS(id, "coefficient"), "Coefficient:",
                            c("Pearson’s coefficient" = "pearson",
                              " Spearman’s rank coefficient" = "spearman"))
         ),
         column(2, 
                # selectInput(NS(id, "by"), "By:",
                #             c("None" = 1,
                #               "Group" = 2,
                #               "Cluster" = 3), selected=1)
                materialSwitch(inputId = NS(id, "clustered"), label = "Show by cluster: ", value = FALSE, status = "primary")
         ),
         column(4, 
                sliderInput(NS(id, "scale"),
                            "Scale:",
                            min = -1,
                            max = 1,
                            step = 0.05,
                            value = c(-1, 1))
         ),
         column(2, 
                downloadButton(NS(id, "download_corr"), "Download Data")
         )
       ),
       add_busy_spinner(spin = "fading-circle", color = "#000000"),
       div(id = "plot-container",
           plotOutput(NS(id, "correlation"), height = "85%")
       )
       
      )
    )
  )
}

