library(shiny)
library(ggplot2)
library(plotly)
library(DT)
library(logger)
library(dplyr)

individualGeneTabUI <- function(id) {
  ns <- NS(id)

  fluidPage(
    titlePanel("Individual Gene Exploration"),
    add_busy_spinner(spin = "fading-circle", color = "#000000"),
    
    fluidRow(
      column(3,
        wellPanel(
          tags$h3("Select Groups & Samples"),
          tags$p("Select which experimental groups and samples you want to analyze."),
          tags$p(tags$b("Note:"), "Your selection here is independent from the main tab."),
          tags$p("You can select any groups and samples for analysis regardless of what you selected in the main tab."),
          uiOutput(ns("toggle_groups_ui")), 
          uiOutput(ns("sample_selection_ui")),
          actionButton(ns("apply_group_selection"), "Apply Group/Sample Selection", 
                       class = "btn-primary", style = "margin-top: 15px;"),
          tags$hr(),
          tags$h3("Select Genes"),
          tags$p("Then, select genes to explore:"),
          selectizeInput(ns("select_gene"), "Select Genes:", 
                      choices = NULL, 
                      multiple = TRUE,
                      options = list(
                        placeholder = "Please select groups & samples first",
                        maxItems = 10,
                        maxOptions = 10000,
                        selectOnTab = TRUE,
                        closeAfterSelect = TRUE
                      ))
        )
      ),
      column(9,
        uiOutput(ns("gene_plots"))
      )
    )
  )
}

individualGeneTabServer <- function(id, dataset, main_dataset) {
  add_group_traces <- function(p, data, ordered_groups, color_mapping) {
    for (grp in ordered_groups) {
      grp_data <- subset(data, Group == grp)
      grp_hover_text <- paste("Replicate: ", grp_data$Sample)

      # Add box plot trace
      p <- add_trace(p,
        data = grp_data, x = ~Group, y = ~Count, type = "box",
        color = I(color_mapping[grp]), name = grp,
        text = grp_hover_text, hoverinfo = "text+y"
      )

      # Add scatter plot points
      grp_colors <- rep(color_mapping[grp], times = nrow(grp_data))

      p <- add_trace(p,
        data = grp_data, x = ~Group, y = ~Count,
        type = "scatter", mode = "markers",
        marker = list(color = I(grp_colors), size = 10, opacity = 0.6),
        text = grp_hover_text, hoverinfo = "text+y",
        showlegend = FALSE
      )
    }
    return(p)
  }

  plot_the_data <- function(data, name, cpm_fpkm_label, gene_name) {
    data$Group <- as.factor(data$Group) # convert to factor for plotting purposes
    ordered_groups <- unique(data$Group)
    upper_limit <- max(data$Count) + 20

    color_mapping <- setNames(object = unique(data$Color), nm = unique(data$Group))
    hover_text <- paste("Replicate: ", data$Sample) # assuming 'SampleName' is the column with individual sample names

    gene_name <- paste(toupper(substring(gene_name, 1, 1)), substring(gene_name, 2), sep = "")

    # Initialize an empty plotly object
    p <- plot_ly()

    p <- add_group_traces(p, data, ordered_groups, color_mapping) %>%
      layout(
        title = list(
          text = paste0("<span style='color:blue;'>", gene_name, "</span>"),
          font = list(size = 20),
          x = 0, # Align title to the left
          xref = "paper"
        ),
        annotations = list(
          list(
            text = name,
            font = list(size = 15, family = "Arial Italic"),
            x = 0.5,
            y = -0.5,
            xref = "paper",
            yref = "paper",
            showarrow = FALSE,
            xanchor = "center",
            yanchor = "bottom"
          )
        ),
        xaxis = list(title = "", categoryorder = "array", categoryarray = ordered_groups),
        yaxis = list(
          title = cpm_fpkm_label,
          zeroline = TRUE,
          zerolinewidth = 2,
          zerolinecolor = "#000000",
          range = c(0, upper_limit)
        ),
        margin = list(t = 100, b = 300, l = 50) # Adjust bottom margin to give space for subtitle
      )

    return(p)
  }

  moduleServer(id, function(input, output, session) {
    # Function to make safe IDs
    make_safe_id <- function(name) {
      gsub("[^A-Za-z0-9_]", "_", name)
    }
    
    # Reactive values for group and sample selection
    all_groups_selected <- reactiveVal(FALSE)
    selected_groups_local <- reactiveVal(c())
    selected_samples_local <- reactiveVal(c())
    selection_applied <- reactiveVal(FALSE)
    
    # Track which samples belong to each group
    samples_by_group <- reactiveVal(list())
    
    # Get complete data from main_dataset instead of filtered data
    complete_data <- reactive({
      req(main_dataset)
      main_dataset$expression_data()
    })

    # Get groups from maintab
    groups <- reactive({
      req(main_dataset)
      main_dataset$groups_data()
    })
    
    # Initialize samples by group mapping
    observe({
      req(groups())
      grp_data <- groups()
      
      # Create a list of samples for each group
      grp_list <- list()
      for (group in unique(grp_data$Group)) {
        samples_in_group <- grp_data %>%
          dplyr::filter(Group == group) %>%
          pull(Sample)
        grp_list[[group]] <- samples_in_group
      }
      
      samples_by_group(grp_list)
    })
    
    # Observe group checkbox changes
    observeEvent(samples_by_group(), {
      # Initialize with all available groups
      selected_groups_local(names(samples_by_group()))
      
      # Initialize with all available samples
      all_samples <- unlist(samples_by_group())
      selected_samples_local(all_samples)
    }, once = TRUE)
    
    # Observe group checkbox changes to update sample selections
    observeEvent(input$toggle_all_groups, {
      req(groups())
      grp_data <- groups()
      all_grps <- unique(grp_data$Group)
      new_state <- !all_groups_selected()
      
      # Update group checkboxes
      for (group in all_grps) {
        group_id <- make_safe_id(paste0("group_", group))
        updateCheckboxInput(session, group_id, value = new_state)
      }
      all_groups_selected(new_state)
      
      # Update sample checkboxes based on group selections
      updateSampleSelections()
    })
    
    # Function to update sample selections based on group selections
    updateSampleSelections <- function() {
      req(samples_by_group())
      grp_data <- groups()
      all_grps <- unique(grp_data$Group)
      grp_list <- samples_by_group()
      
      # Get currently selected groups
      selected_grps <- Filter(function(group) 
        isTRUE(input[[make_safe_id(paste0("group_", group))]]), all_grps)
      
      # For each group, update its samples
      for (group in all_grps) {
        group_selected <- group %in% selected_grps
        samples_in_group <- grp_list[[group]]
        
        if (!is.null(samples_in_group)) {
          for (sample in samples_in_group) {
            sample_id <- make_safe_id(paste0("sample_", sample))
            updateCheckboxInput(session, sample_id, value = group_selected)
          }
        }
      }
    }
    
    # When a group checkbox changes, update its samples
    observeEvent(lapply(names(samples_by_group()), function(group) {
      input[[make_safe_id(paste0("group_", group))]]
    }), {
      updateSampleSelections()
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
    
    # Dynamically render the group selection UI
    output$toggle_groups_ui <- renderUI({
      req(groups())
      actionButton(session$ns("toggle_all_groups"), 
                  label = ifelse(all_groups_selected(), 
                                "Deselect All Groups", 
                                "Select All Groups"))
    })
    
    # Dynamically render the groups and samples
    output$sample_selection_ui <- renderUI({
      req(groups(), samples_by_group())
      grp_data <- groups()
      all_groups <- unique(grp_data$Group)
      grp_list <- samples_by_group()
      
      tagList(
        tags$h4("Toggle Groups"),
        tags$div(style = "max-height: 400px; overflow-y: auto;",
          lapply(all_groups, function(group) {
            group_id <- make_safe_id(paste0("group_", group))
            samples_in_group <- grp_list[[group]]
            
            tagList(
              checkboxInput(session$ns(group_id), group, value = TRUE),
              # Show selectable samples for this group
              conditionalPanel(
                condition = sprintf("input['%s'] == true", session$ns(group_id)),
                tags$div(style = "margin-left: 20px; margin-bottom: 10px;",
                  lapply(samples_in_group, function(sample) {
                    sample_id <- make_safe_id(paste0("sample_", sample))
                    checkboxInput(session$ns(sample_id), sample, value = TRUE)
                  })
                )
              )
            )
          })
        )
      )
    })
    
    # Apply group selection
    observeEvent(input$apply_group_selection, {
      req(groups(), complete_data(), samples_by_group())
      grp_data <- groups()
      all_groups <- unique(grp_data$Group)
      grp_list <- samples_by_group()
      
      # Get selected groups
      selected_groups_list <- Filter(function(group) 
        isTRUE(input[[make_safe_id(paste0("group_", group))]]), all_groups)
      
      if (length(selected_groups_list) == 0) {
        showNotification("No groups selected. Please select at least one group.", type = "error")
        return(NULL)
      }
      
      selected_groups_local(selected_groups_list)
      
      # Get all samples in selected groups
      samples_in_selected_groups <- unlist(grp_list[selected_groups_list])
      
      # Get selected samples from checkboxes
      selected_samples_list <- samples_in_selected_groups %>%
        Filter(function(sample) 
          isTRUE(input[[make_safe_id(paste0("sample_", sample))]]), .)
      
      if (length(selected_samples_list) == 0) {
        showNotification("No samples selected. Please select samples to display.", type = "error")
        return(NULL)
      }
      
      selected_samples_local(selected_samples_list)
      selection_applied(TRUE)
      
      # Show confirmation
      showNotification("Group and sample selection applied", type = "message")
      
      # Update gene dropdown
      updateGeneList()
    })
    
    # Update gene list after selection is applied - this now uses complete_data instead of filtered_data
    updateGeneList <- function() {
      req(complete_data(), selection_applied())
      genes <- complete_data()$Gene_Symbol
      updateSelectizeInput(session, "select_gene", 
                        choices = genes,
                        selected = if(length(genes) > 0) genes[1] else NULL,
                        server = TRUE)
    }

    # Filter the data based on the selected gene and custom group/sample selection
    # This now uses complete_data instead of filtered_data from main tab
    filtered_data <- reactive({
      req(input$select_gene, complete_data(), selected_groups_local(), selected_samples_local(), selection_applied())
      
      # Selected groups and samples
      selected_grps <- selected_groups_local()
      selected_samples <- selected_samples_local()
      
      # Get available columns in the data
      available_columns <- colnames(complete_data())
      
      # Filter to only include samples that exist in the data
      valid_samples <- intersect(selected_samples, available_columns)
      
      # Filter the data based on the selected gene
      res <- complete_data() %>%
        filter(Gene_Symbol %in% input$select_gene)
      
      # Select only Symbol, Gene_Symbol and valid samples 
      if (length(valid_samples) > 0) {
        res <- res %>% select(Symbol, Gene_Symbol, all_of(valid_samples))
      } else {
        res <- res %>% select(Symbol, Gene_Symbol)
      }

      # Reshape the data
      res_long <- res %>%
        pivot_longer(cols = -c(Symbol, Gene_Symbol), names_to = "Sample", values_to = "Count")

      # Merge the data with groups, filtering for selected groups
      merged_data <- res_long %>%
        left_join(groups(), by = "Sample") %>%
        filter(Group %in% selected_grps)

      return(merged_data)
    })

    # Render the gene plots
    output$gene_plots <- renderUI({
      if (!selection_applied()) {
        return(tags$div(
          class = "alert alert-info",
          tags$h4("Please Select Groups and Samples"),
          tags$p("Select groups and samples on the left panel, then click 'Apply Group/Sample Selection' to load available genes.")
        ))
      }
      
      if (is.null(input$select_gene) || length(input$select_gene) == 0) {
        return(tags$div(
          class = "alert alert-warning",
          tags$h4("No Genes Selected"),
          tags$p("Please select one or more genes from the dropdown menu on the left.")
        ))
      }
      
      req(filtered_data())
      ns <- session$ns
      gene_names <- unique(filtered_data()$Gene_Symbol)
      
      if (length(gene_names) == 0) {
        return(tags$div(
          class = "alert alert-warning",
          tags$h4("No Data for Selected Genes"),
          tags$p("There is no data for the selected genes in the chosen groups and samples.")
        ))
      }
      
      plot_output_list <- lapply(gene_names, function(gene) {
        tags$div(
          class = "well",
          style = "margin-bottom: 20px;",
          plotlyOutput(ns(paste0("plot_", gene)), height = "500px", width = "100%")
        )
      })
      do.call(tagList, plot_output_list)
    })

    observe({
      req(filtered_data())
      gene_data <- filtered_data()
      gene_names <- unique(gene_data$Gene_Symbol)

      lapply(gene_names, function(gene) {
        output[[paste0("plot_", gene)]] <- renderPlotly({
          gene_specific_data <- gene_data %>% filter(Gene_Symbol == gene)
          p <- plot_the_data(gene_specific_data,
            name = "Gene Expression Data",
            cpm_fpkm_label = "Expression Level", gene_name = gene
          )
          p %>% layout(
            autosize = TRUE,
            margin = list(l = 50, r = 50, t = 80, b = 80)
          )
        })
      })
    })
  })
}
