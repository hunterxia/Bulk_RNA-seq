library(ggplot2)
library(corrr)
PCACorrelationTabServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    data <- reactive({
      req(input$expressions)
      ext <- tools::file_ext(input$expressions$name)
      switch(ext,
             csv = vroom::vroom(input$expressions$datapath, delim = ","),
             tsv = vroom::vroom(input$expressions$datapath, delim = "\t"),
             validate("Invalid file; Please upload a .csv or .tsv file")
      )
    })
    
    output$pca <- renderPlot({
      req(input$expressions)
      df <- data()
      colnames <- colnames(df)
      df_transposed <- t(df[,3:ncol(df)])
      pca_comp <- prcomp(df_transposed, scale. = TRUE, center = TRUE)
      percentVar <- pca_comp$sdev^2/sum(pca_comp$sdev^2)
      pca_df <- data.frame(PC1= pca_comp$x[,1], PC2= pca_comp$x[,2], sampleLab = rownames(pca_comp$x), check.names = FALSE)
      
      pca_df %>%
        ggplot(aes(x = PC1,y = PC2, label = sampleLab)) +
        geom_point(size=5)
    }, res = 96)
    
    
    output$correlation <- renderPlot({
      req(input$expressions)
      df <- data()
      numerical_data <- df[,3:ncol(df)]
      data_normalized <- scale(numerical_data)
      corr_matrix <- cor(data_normalized, method = input$coefficient)
      ggcorrplot(corr_matrix)
    }, res = 96)
    
    
  })
}

PCACorrelationTabUI <- function(id) {
  fluidPage(
    titlePanel("PCA and Correlation"),
    fluidRow(
      column(12, 
        fileInput(NS(id, "expressions"), NULL, accept = c(".csv", ".tsv"))
      )
    ),
    fluidRow(
      column(6, 
       selectInput(NS(id, "pc"), "Which PCs:",
                   c("PC1 vs. PC2" = "pearson")),
        plotOutput(NS(id, "pca"))
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

