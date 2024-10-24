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
      
      numericInput(ns("expression_cutoff"), "Expression Level Cutoff", value = 0),
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
    expression_data <- reactiveVal()       # Stores the uploaded expression data
    groups_data <- reactiveVal()           # Stores the uploaded groups data
    filtered_data <- reactiveVal()         # Stores the filtered expression data
    selected_groups <- reactiveVal()       # Stores the selected groups
    selected_samples_data <- reactiveVal() # Stores the selected samples
    
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
      groups_data(grp_data)  # Store group data in reactive variable
      
      # Dynamically render the groups and samples with indentation
      output$sample_selection_ui <- renderUI({
        req(groups_data())
        grp_data <- groups_data()
        
        tagList(
          tags$h4("Toggle Groups"),  # Add a title for the toggle groups section
          
          lapply(unique(grp_data[['Group']]), function(group) {
            samples_in_group <- grp_data %>%
              dplyr::filter(Group == group) %>%
              pull(Sample)
            
            tagList(
              checkboxInput(ns(paste0("group_", group)), group, value = FALSE),  # Group checkbox
              conditionalPanel(
                condition = sprintf("input['%s'] == true", ns(paste0("group_", group))),
                tagList(
                  lapply(samples_in_group, function(sample) {
                    div(style = "margin-left: 20px;",  # Indent the samples
                        checkboxInput(ns(paste0("sample_", sample)), sample, value = TRUE)
                    )
                  })
                )
              )
            )
          })
        )
      })
    })
    
    # Reading and storing expression data
    observeEvent(input$expression_input, {
      req(input$expression_input)
      data <- read_data(input$expression_input)
      # Convert the numeric columns (third to last) to numeric types
      data[, 3:ncol(data)] <- lapply(data[, 3:ncol(data)], as.numeric)
      names(data) <- gsub(" ", "_", names(data))  # Clean column names
      expression_data(data)  # Store expression data in reactive variable
    })
    
    # Apply filters to expression data based on selected samples and cutoff values
    observeEvent(input$apply_filters, {
      req(expression_data(), input$expression_cutoff)
      expr_data <- expression_data()
      
      # Collect selected samples from the checkbox inputs
      selected_samples <- groups_data() %>%
        pull(Sample) %>%
        Filter(function(sample) input[[paste0("sample_", sample)]], .)
      
      if (length(selected_samples) == 0) {
        showNotification("No samples selected. Please select samples to apply filters.", type = "error")
        return(NULL)
      }
      
      # Update selected samples data reactive value
      selected_samples_data(selected_samples)
      cutoff <- input$expression_cutoff
      
      # Filtering the expression data based on selected samples and cutoff
      filtered_expr_data <- expr_data %>%
        select(c('Symbol', 'Gene_Symbol', all_of(selected_samples))) %>%
        dplyr::filter(rowSums(.[selected_samples] >= cutoff) >= input$sample_count_cutoff)
      
      if (nrow(filtered_expr_data) == 0) {
        showNotification("No data matches the filtering criteria. Please adjust your filters.", type = "error")
        return(NULL)
      }
      
      filtered_data(filtered_expr_data)  # Store filtered data in reactive variable
    })
    
    # Render the filtered data table
    output$resultsTable <- renderDT({
      req(filtered_data())
      filtered_data()
    }, options = list(pageLength = 25))
    
    # Download handler for the filtered data
    output$download_filtered_data <- downloadHandler(
      filename = function() { "filtered_data.csv" },
      content = function(file) {
        req(filtered_data())
        write.csv(filtered_data(), file, row.names = FALSE)
      }
    )
    
    # Return a list of reactive values for other components to access
    return(list(
      expression_data = expression_data,
      groups_data = groups_data,
      filtered_data = filtered_data,
      selected_groups = selected_groups,
      selected_samples_data = selected_samples_data
    ))
  })
}
