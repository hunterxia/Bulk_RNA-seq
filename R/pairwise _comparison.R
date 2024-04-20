library(shiny)
library(ggplot2)
library(DT)

PairwiseComparisonTabUI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    titlePanel("Pairwise Comparison of Gene Expression"),
    sidebarLayout(
      sidebarPanel(
        selectInput(ns("select_group1"), "Select Group 1", choices = NULL ),
        selectInput(ns("select_group2"), "Select Group 2", choices = NULL),
        selectInput(ns("test_type"), "Select Test Type",
                    choices = c("Paired t-test" = "paired",
                                "Unpaired t-test" = "unpaired",
                                "Non-parametric test" = "nonparametric")),
        numericInput(ns("fc_cutoff"), "Fold Change Cutoff", value = 1),
        numericInput(ns("pvalue_cutoff"), "P-value Cutoff", value = 0.05),
        checkboxInput(ns("adjust_pvalue"), "Use Adjusted P-Value", value = FALSE),
        actionButton(ns("run_analysis"), "Run Analysis"),
        downloadButton(ns("download_all"), "Download Full DEG Output"),
        downloadButton(ns("download_filtered"), "Download Filtered DEGs")
      ),
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Volcano Plot", plotOutput(ns("volcano_plot"))),
                    tabPanel("MA Plot", plotOutput(ns("ma_plot"))),
                    tabPanel("DEG Table", DTOutput(ns("deg_table")))
        )
      )
    )
  )
}



PairwiseComparisonTabServer <- function(id,dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Placeholder reactive value for storing expression data and results
    expression_data <- reactiveVal()
    deg_results <- reactiveVal()
    
    observe({
      updateSelectInput(session, ns("select_group1"), choices = names(dataset$groups_data()))
      updateSelectInput(session, ns("select_group2"), choices = names(dataset$groups_data()))
    })
    
    observeEvent(input$run_analysis, {

      if(input$test_type == "paired") {
        # Run paired t-test
      } else if(input$test_type == "unpaired") {
        # Run unpaired t-test
      } else {
        # Run non-parametric test
      }
      
      # Store DEGs, calculate Log2FC, p-values, etc.
      # deg_results() <- some calculation
      
      
      output$volcano_plot <- renderPlot({
        # plot volcano plot using deg_results()
      })
      
      output$ma_plot <- renderPlot({
        # plot MA plot using deg_results()
      })
      
      output$deg_table <- renderDT({
        # create a DataTable from deg_results()
      })
      
      # Download handlers
      output$download_all <- downloadHandler(
        filename = function() { "full_deg_output.csv" },
        content = function(file) {
          # write.csv(deg_results(), file)
        }
      )
      
      output$download_filtered <- downloadHandler(
        filename = function() { "filtered_deg_output.csv" },
        content = function(file) {
        }
      )
    })
  })
}