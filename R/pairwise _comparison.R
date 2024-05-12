library(shiny)
library(ggplot2)
library(DT)
library(dplyr)

PairwiseComparisonTabUI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    titlePanel("Pairwise Comparison of Gene Expression"),
    sidebarLayout(
      sidebarPanel(
        selectInput(ns("select_group1"), "Select Group 1", choices = NULL),
        selectInput(ns("select_group2"), "Select Group 2", choices = NULL),
        selectInput(ns("test_type"), "Select Test Type",
                    choices = c("Paired t-test" = "paired",
                                "Unpaired t-test" = "unpaired")),
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


perform_DEG_analysis <- function(exp_data, group1_samples, group2_samples, test_type, fc_cutoff, pvalue_cutoff, adjust_pvalue) {
  # Filter expression data to include only the relevant samples
  group1_data <- exp_data[, colnames(exp_data) %in% group1_samples, drop = FALSE]
  group2_data <- exp_data[, colnames(exp_data) %in% group2_samples, drop = FALSE]
  
  # Ensure the data is numeric
  group1_data <- data.frame(lapply(group1_data, as.numeric))
  group2_data <- data.frame(lapply(group2_data, as.numeric))
  
  # Calculate log2 fold changes
  log2_fold_changes <- log2(rowMeans(group2_data + 1, na.rm = TRUE)) - log2(rowMeans(group1_data + 1, na.rm = TRUE))
  
  # Calculate p-values for each gene
  p_values <- vector("numeric", length = nrow(exp_data))
  for (i in 1:nrow(exp_data)) {
    gene_data1 <- group1_data[i, , drop = FALSE]  # Ensure it remains as a dataframe
    gene_data2 <- group2_data[i, , drop = FALSE]  # Ensure it remains as a dataframe
    
    if (test_type == "paired") {
      p_values[i] <- t.test(as.numeric(gene_data1), as.numeric(gene_data2), paired = TRUE)$p.value
    } else {
      p_values[i] <- t.test(as.numeric(gene_data1), as.numeric(gene_data2), paired = FALSE)$p.value
    }
  }
  
  results <- data.frame(
    Gene = rownames(exp_data),
    Log2FoldChange = log2_fold_changes,
    PValue = p_values
  )
  
  # Adjust p-values if requested
  if (adjust_pvalue) {
    results$AdjPValue <- p.adjust(results$PValue, method = "BH")
  } else {
    results$AdjPValue <- results$PValue
  }
  
  # Filter results based on cutoffs
  results <- results[abs(results$Log2FoldChange) >= log2(fc_cutoff) & results$AdjPValue <= pvalue_cutoff, ]
  return(results)
}


# Volcano plot generation function
generate_volcano_plot <- function(results) {
  ggplot(results, aes(x = Log2FoldChange, y = -log10(PValue), color = AdjPValue < 0.05)) +
    geom_point(alpha = 0.8) +
    labs(x = "Log2 Fold Change", y = "-log10(P-value)", title = "Volcano Plot") +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal()
}

# MA plot generation function
generate_ma_plot <- function(exp_data, group1_samples, group2_samples) {
  A <- (log2(rowMeans(exp_data[, group1_samples, drop = FALSE] + 1)) + log2(rowMeans(exp_data[, group2_samples, drop = FALSE] + 1))) / 2
  M <- log2(rowMeans(exp_data[, group2_samples, drop = FALSE] + 1)) - log2(rowMeans(exp_data[, group1_samples, drop = FALSE] + 1))
  ggplot(data.frame(A = A, M = M), aes(x = A, y = M)) +
    geom_point(alpha = 0.5) +
    labs(x = "Average Log Intensity (A)", y = "Log Fold Change (M)", title = "MA Plot") +
    theme_minimal()
}

PairwiseComparisonTabServer <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Access the processed data from the dataset reactive
    expression_data <- reactive(dataset$expression_data())
    groups_data <- reactive(dataset$groups_data())
    
    group_names_vector <- reactive({
      unique(groups_data()$Exp_Grp)
    })
    
    observe({
      grp_names <- group_names_vector()
      if (is.null(grp_names)) return()
      updateSelectInput(session, "select_group1", choices = grp_names, selected = grp_names[1])
      updateSelectInput(session, "select_group2", choices = grp_names, selected = grp_names[length(grp_names)])
    })
    
    observeEvent(input$run_analysis, {
      req(expression_data(), groups_data())
      exp_data <- expression_data()
      grp_data <- groups_data()
      
      # Map group selections to sample names
      group1_samples <- grp_data$HQ_samples[grp_data$Exp_Grp == input$select_group1]
      group2_samples <- grp_data$HQ_samples[grp_data$Exp_Grp == input$select_group2]
      
      
      # Perform the DEG analysis
      results <- perform_DEG_analysis(exp_data, group1_samples, group2_samples,
                                      input$test_type, input$fc_cutoff, input$pvalue_cutoff, input$adjust_pvalue)
      
      # Update outputs: Volcano Plot, MA Plot, DEG Table
      output$volcano_plot <- renderPlot({ generate_volcano_plot(results) })
      output$ma_plot <- renderPlot({ generate_ma_plot(exp_data, group1_samples, group2_samples) })
      output$deg_table <- renderDT({ datatable(results, options = list(pageLength = 10)) })
    })
    
    # Download handlers for DEG results
    output$download_all <- downloadHandler(
      filename = function() { "full_deg_output.csv" },
      content = function(file) {
        write.csv(results(), file, row.names = FALSE)
      }
    )
    output$download_filtered <- downloadHandler(
      filename = function() { "filtered_deg_output.csv" },
      content = function(file) {
        write.csv(subset(results(), AdjPValue <= input$pvalue_cutoff), file, row.names = FALSE)
      }
    )
  })
}
