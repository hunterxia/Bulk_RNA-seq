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
    expression_data <- reactiveVal()
    groups_data <- reactiveVal()
    filtered_data <- reactiveVal()
    selected_groups <- reactiveVal()
    selected_samples_data <- reactiveVal()
    
    # Function to make safe IDs
    make_safe_id <- function(name) {
      gsub("[^A-Za-z0-9_]", "_", name)
    }
    
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
      
      # Dynamically render the groups and samples with sanitized IDs
      output$sample_selection_ui <- renderUI({
        req(groups_data())
        grp_data <- groups_data()
        
        tagList(
          tags$h4("Toggle Groups"),
          
          lapply(unique(grp_data[['Group']]), function(group) {
            group_id <- make_safe_id(paste0("group_", group))
            samples_in_group <- grp_data %>%
              dplyr::filter(Group == group) %>%
              pull(Sample)
            
            tagList(
              checkboxInput(ns(group_id), group, value = FALSE),  # Group checkbox
              conditionalPanel(
                condition = sprintf("input['%s'] == true", ns(group_id)),
                tagList(
                  lapply(samples_in_group, function(sample) {
                    sample_id <- make_safe_id(paste0("sample_", sample))
                    div(style = "margin-left: 20px;",
                        checkboxInput(ns(sample_id), sample, value = TRUE)
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
      data[, 3:ncol(data)] <- lapply(data[, 3:ncol(data)], as.numeric)
      names(data) <- gsub(" ", "_", names(data))
      expression_data(data)
    })
    
    # Apply filters to expression data
    observeEvent(input$apply_filters, {
      req(expression_data(), input$expression_cutoff)
      expr_data <- expression_data()
      grp_data <- groups_data()
      
      # Get all groups
      all_groups <- unique(grp_data[['Group']])
      
      # Get selected groups
      selected_groups_list <- Filter(function(group) isTRUE(input[[make_safe_id(paste0("group_", group))]]), all_groups)
      selected_groups(selected_groups_list)
      
      if (length(selected_groups_list) == 0) {
        showNotification("No groups selected. Please select at least one group.", type = "error")
        return(NULL)
      }
      
      # Get samples in selected groups
      samples_in_selected_groups <- grp_data %>%
        filter(Group %in% selected_groups_list) %>%
        pull(Sample)
      
      # Get selected samples
      selected_samples_list <- samples_in_selected_groups %>%
        Filter(function(sample) isTRUE(input[[make_safe_id(paste0("sample_", sample))]]), .)
      
      if (length(selected_samples_list) == 0) {
        showNotification("No samples selected. Please select samples to apply filters.", type = "error")
        return(NULL)
      }
      
      selected_samples_data(selected_samples_list)
      cutoff <- input$expression_cutoff
      
      # Filtering the expression data
      filtered_expr_data <- expr_data %>%
        select(c('Symbol', 'Gene_Symbol', all_of(selected_samples_list))) %>%
        dplyr::filter(rowSums(.[selected_samples_list] >= cutoff) >= input$sample_count_cutoff)
      
      if (nrow(filtered_expr_data) == 0) {
        showNotification("No data matches the filtering criteria. Please adjust your filters.", type = "error")
        return(NULL)
      }
      
      filtered_data(filtered_expr_data)
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
    
    # Return reactive values
    return(list(
      expression_data = expression_data,
      groups_data = groups_data,
      filtered_data = filtered_data,
      selected_groups = selected_groups,
      selected_samples_data = selected_samples_data
    ))
  })
}
