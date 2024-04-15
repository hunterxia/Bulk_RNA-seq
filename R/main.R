library(shiny)
library(DT)
library(readxl)
library(readr)
library(dplyr)

mainTabUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      fileInput(ns("expression_input"), "Upload Expression Table", 
                accept = c(".xlsx", ".csv")),
      fileInput(ns("groups_input"), "Upload Groups Table", 
                accept = c(".xlsx", ".csv")),
      uiOutput(ns("toggle_groups_ui")),
      numericInput(ns("expression_cutoff"), "Expression Cutoff", min = 0, max = 100, value = 1, step = 0.1),
      uiOutput(ns("sample_selection_ui")),
      actionButton(ns("apply_filters"), "Apply Filters")
      # Sidebar content ends here
    ),
    mainPanel(
      DTOutput(ns("resultsTable"))
      # Main panel content ends here
    )
  )
}



mainTabServer <- function(id) {
    moduleServer(id, function(input, output, session) {
      ns <- session$ns
      
      # Reactive values to store data
      expression_data <- reactiveVal()
      groups_data <- reactiveVal()
      filtered_data <- reactiveVal()  # For storing filtered data
      
      # Define a function to read data based on file type
      read_data <- function(file) {
        if (grepl("\\.xlsx$", file$name)) {
          read_excel(file$datapath)
        } else if (grepl("\\.csv$", file$name)) {
          read_csv(file$datapath, col_types = cols(.default = "c"))
        } else {
          stop("Unsupported file type")
        }
      }
      
      # Reading and storing group data
      observeEvent(input$groups_input, {
        req(input$groups_input)
        grp_data <- read_data(input$groups_input)
        groups_data(grp_data)
        
        output$toggle_groups_ui <- renderUI({
          checkboxGroupInput(ns("toggle_groups"), "Toggle Groups", choices = unique(grp_data[['Exp_Grp']]))
        })
        output$sample_selection_ui <- renderUI({
          checkboxGroupInput(ns("sample_selection"), "Select Samples", choices = grp_data$HQ_samples)
        })
      })
      
      # Reading and storing expression data
      observeEvent(input$expression_input, {
        req(input$expression_input)
        data <- read_data(input$expression_input)
        names(data) <- gsub(" ", "_", names(data))
        expression_data(data)
      })
      
      # Apply filters based on user inputs
      observeEvent(input$apply_filters, {
        req(expression_data(), groups_data(), input$sample_selection, input$expression_cutoff)
        
        expr_data <- expression_data()
        selected_samples <- input$sample_selection
        selected_groups <- input$toggle_groups
        cutoff <- input$expression_cutoff
        
        valid_samples <- intersect(selected_samples, groups_data() %>% filter(Exp_Grp %in% selected_groups) %>% .$HQ_samples)
        
        filtered_expr_data <- expr_data %>%
          select(c('Symbols', 'Genes', all_of(valid_samples))) %>%
          filter(rowSums(.[valid_samples] >= cutoff) >= 1)
        
        filtered_data(filtered_expr_data)
      })
      
      # Render the filtered data table
      output$resultsTable <- renderDT({
        groups_data()
      }, options = list(pageLength = 25))
      
      # return a list of data for other modules to access
      return(list(
        expression_data = expression_data,
        groups_data = groups_data,
        filtered_data = filtered_data
      ))
    })
  }
  