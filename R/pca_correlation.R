library(ggplot2)
library(corrr)
library(plotly)
library(ggcorrplot)

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
      data_normalized <- scale(numerical_data)
      corr_matrix <- cor(data_normalized, method = input$coefficient)
      ggcorrplot(corr_matrix)
    }
    
    
    output$pca <- renderPlotly({
      req(dataset$expression_data())
      create_pca_plot()
    })
    
    output$correlation <- renderPlot({
      req(dataset$expression_data())
      create_correlation_plot()
    }, res = 96)
    
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
         )
       ),
       plotlyOutput(NS(id, "pca"), width = "100%", height = "700px")
      ),
      column(6, 
       selectInput(NS(id, "coefficient"), "Coefficient:",
                   c("Pearson’s coefficient" = "pearson",
                     " Spearman’s rank coefficient" = "spearman")),
        plotOutput(NS(id, "correlation"))
      )
    )
  )
}

