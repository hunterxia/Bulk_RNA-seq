library(shiny)
library(ggplot2)
library(plotly)
library(DT)
library(logger)

individualGeneTabUI <- function(id) {
  ns <- NS(id)

  fluidPage(
    titlePanel("Individual Gene"),
    add_busy_spinner(spin = "fading-circle", color = "#000000"),
    sidebarLayout(
      sidebarPanel(
        selectInput(ns("select_gene"), "Select Gene", choices = c("Please wait..." = ""), multiple = TRUE)
      ),
      mainPanel(
        uiOutput(ns("gene_plots")),
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
    # get genes from the clustering data
    clustering_data <- reactive({
      req(dataset)
      dataset$selected_variable_genes()
    })

    # get data from maintab
    main_filtered_data <- reactive({
      req(main_dataset)
      main_dataset$filtered_data()
    })

    # get groups from maintab
    groups <- reactive({
      req(main_dataset)
      main_dataset$groups_data()
    })

    # Update the select gene input
    observe({
      req(main_filtered_data())
      genes <- main_filtered_data()$Gene_Symbol
      updateSelectInput(session, "select_gene", choices = genes, selected = genes[1])
    })

    # Filter the data based on the selected gene
    filtered_data <- reactive({
      req(input$select_gene, main_filtered_data())

      # Filter the data based on the selected gene
      res <- main_filtered_data() %>%
        filter(Gene_Symbol %in% input$select_gene)


      # Reshape the data
      res_long <- res %>%
        pivot_longer(cols = -c(Symbol, Gene_Symbol), names_to = "Sample", values_to = "Count")

      # Merge the data with groups
      merged_data <- res_long %>%
        left_join(groups(), by = "Sample")

      return(merged_data)
    })

    # Render the gene plots
    output$gene_plots <- renderUI({
      req(filtered_data())
      ns <- session$ns
      gene_names <- unique(filtered_data()$Gene_Symbol)
      plot_output_list <- lapply(gene_names, function(gene) {
        plotlyOutput(ns(paste0("plot_", gene)), height = "600px", width = "auto")
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
