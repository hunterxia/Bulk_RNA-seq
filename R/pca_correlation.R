library(ggplot2)
library(corrr)
library(plotly)
library(ggcorrplot)
library(reshape2)

PCACorrelationTabServer <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    create_pca_plot <- function(){
      req(input$pc_x, input$pc_y)
      pc_x <- as.numeric(input$pc_x)
      pc_y <- as.numeric(input$pc_y)
      df <- dataset$expression_data()
      colnames <- colnames(df)
      df_transposed <- t(df[,3:ncol(df)])
      pca_comp <- prcomp(df_transposed, scale. = TRUE, center = TRUE)
      percentVar <- pca_comp$sdev^2/sum(pca_comp$sdev^2)
      pca_df <- data.frame(pca_comp$x, sampleLab = rownames(pca_comp$x), check.names = FALSE)
      
      p <- pca_df %>%
        ggplot(aes(x = pca_df[, pc_x],y = pca_df[, pc_y], label = sampleLab, color = sampleLab, 
                   text = paste("Sample:", sampleLab, "<br>PC", pc_x ,":", pca_df[, pc_x], "<br>PC", pc_y  ,":", pca_df[, pc_y]))) +
        geom_point(size=5) +
        labs(x=paste0("PC", pc_x, ": ", round(percentVar[pc_x]*100,1), "% variance"),
             y=paste0("PC", pc_y, ": ", round(percentVar[pc_y]*100,1), "% variance")) +
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
    
    create_correlation_plot <- function(){
      df <- dataset$expression_data()
      numerical_data <- df[,3:ncol(df)]
      #data_normalized <- scale(numerical_data)
      corr_matrix <- cor(data_normalized, method = input$coefficient)
      #ggcorrplot(corr_matrix)
      corr_melted <- melt(corr_matrix)
      corr_melted$Var1 <- factor(corr_melted$Var1, levels = rev(unique(corr_melted$Var1)))
      corr_melted$Var2 <- factor(corr_melted$Var2, levels = unique(corr_melted$Var2))
      
      custom_min <- input$scale[1]
      custom_max <- input$scale[2]
      colors <- colorRampPalette(c("blue", "white", "red"))(5)
      breaks <- round(seq(custom_min, custom_max, length.out = 5), 2)
      
      # Plot using ggplot2
      ggplot(data = corr_melted, aes(Var1, Var2, fill = value)) +
        geom_tile() + # Create a heatmap
        #geom_text(aes(label = round(value, 2)), size = 3) + # Add correlation coefficients
        scale_fill_gradientn(colors = colors,
                             values = scales::rescale(breaks),
                             limits = c(custom_min, custom_max),
                             oob = scales::squish,
                             name = "Correlation") +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 5),
              axis.text.y = element_text(size = 5), 
              axis.title = element_blank())
    }
    
    
    output$pca <- renderPlotly({
      req(dataset$expression_data())
      create_pca_plot()
    })
    
    output$correlation <- renderPlotly({
      req(dataset$expression_data())
      create_correlation_plot()
    })
    
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
    observeEvent(dataset$expression_data(), {
      df <- dataset$expression_data()
      df_transposed <- t(df[,3:ncol(df)])
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
         column(3, 
                downloadButton(NS(id, "download_pca"), "Download Data")
         )
       ),
       plotlyOutput(NS(id, "pca"), width = "100%", height = "700px")
      ),
      column(6, 
       fluidRow(
         column(3, 
                selectInput(NS(id, "coefficient"), "Coefficient:",
                            c("Pearson’s coefficient" = "pearson",
                              " Spearman’s rank coefficient" = "spearman"))
         ),
         column(4, 
                sliderInput(NS(id, "scale"),
                            "Scale:",
                            min = -1,
                            max = 1,
                            step = 0.05,
                            value = c(-1, 1))
         ),
         column(3, 
                downloadButton(NS(id, "download_corr"), "Download Data")
         )
       ),
       plotlyOutput(NS(id, "correlation"), height = "100%")
      )
    )
  )
}

