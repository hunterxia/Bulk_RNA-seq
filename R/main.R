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
      uiOutput(ns("sample_selection_ui")),  
      
      numericInput(ns("expression_cutoff"), "Expression Level Cutoff", value = 1),
      numericInput(ns("sample_count_cutoff"), "Minimum Sample Count", value = 1, min = 1),
      
      actionButton(ns("apply_filters"), "Apply Filters"),
      
      downloadButton(ns("download_filtered_data"), "Download Filtered Data")
    ),
    mainPanel(
      DTOutput(ns("resultsTable"))
    )
  )
}



mainTabServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values to store data
    expression_data <- reactiveVal()
    groups_data <- reactiveVal()
    filtered_data <- reactiveVal()  
    selected_groups <- reactiveVal()
    selected_samples_data <- reactiveVal()
    
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
        checkboxGroupInput(ns("toggle_groups"), "Toggle Groups", choices = unique(grp_data[['Exp_Grp']]), selected = NULL)
      })
    })
    
    # Dynamic sample selection based on selected groups
    observe({
      req(input$toggle_groups)
      selected_groups(input$toggle_groups)
    })
    
    observe({
      req(groups_data())
      grp_data <- groups_data()
      selected_groups <- input$toggle_groups
      
      
      if (length(selected_groups) > 0) {
        samples_in_selected_groups <- grp_data %>%
          filter(Exp_Grp %in% selected_groups) %>%
          pull(HQ_samples) %>%
          unique()
        
        output$sample_selection_ui <- renderUI({
          checkboxGroupInput(ns("sample_selection"), "Select Samples", choices = samples_in_selected_groups, selected = samples_in_selected_groups)
        })
      } else {
        
        output$sample_selection_ui <- renderUI({
          
        })
      }
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
      req(expression_data(), input$sample_selection, input$expression_cutoff)
      
      expr_data <- expression_data()
      selected_samples <- input$sample_selection
      # update reactive variable
      selected_samples_data(selected_samples)
      cutoff <- input$expression_cutoff
      
      # Filtering the expression data
      filtered_expr_data <- expr_data %>%
        select(c('Symbols', 'Genes', all_of(selected_samples))) %>%
        filter(rowSums(.[selected_samples] >= cutoff) >= input$sample_count_cutoff)
      
      filtered_data(filtered_expr_data)
    })
    
    # Render the filtered data table
    output$resultsTable <- renderDT({
      filtered_data()
    }, options = list(pageLength = 25))
    
    output$download_filtered_data <- downloadHandler(
      filename = function() { "filtered_data.csv" },
      content = function(file) {
        req(filtered_data())
        write.csv(filtered_data(), file, row.names = FALSE)
      }
    )
    
    # Return a list of data for other modules to access
    return(list(
      expression_data = expression_data,
      groups_data = groups_data,
      filtered_data = filtered_data,
      selected_groups = selected_groups,
      selected_samples_data = selected_samples_data
    ))
  })
}