library(shiny)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)
library(DESeq2)

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
                                "Unpaired t-test" = "unpaired",
                                "Deseq" = "deseq")),
        numericInput(ns("fc_cutoff"), "Fold Change Cutoff", value = 1),
        numericInput(ns("pvalue_cutoff"), "P-value Cutoff", value = 0.05),
        checkboxInput(ns("adjust_pvalue"), "Use Adjusted P-Value", value = FALSE),
        actionButton(ns("run_analysis"), "Run Analysis"),
        downloadButton(ns("download_all"), "Download Full DEG Output"),
        downloadButton(ns("download_filtered"), "Download Filtered DEGs")
      ),
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Volcano Plot",
                             textOutput(ns("volcano_gene_counts")),
                             plotlyOutput(ns("volcano_plot"), width = "700px", height = "700px")),
                    tabPanel("MA Plot",
                             textOutput(ns("ma_gene_counts")),
                             plotlyOutput(ns("ma_plot"), width = "700px", height = "700px")),
                    tabPanel("Positive DEG Table", DTOutput(ns("pos_deg_table"))),
                    tabPanel("Negative DEG Table", DTOutput(ns("neg_deg_table")))
        )
      )
    )
  )
}

perform_DEG_analysis <- function(exp_data, group1_samples, group2_samples, test_type, fc_cutoff, pvalue_cutoff, adjust_pvalue) {
  if (test_type == "deseq") {
    exp_data_filtered <- exp_data[, c(group1_samples, group2_samples)]
    exp_data_filtered <- round(exp_data_filtered)
    
    condition <- factor(c(rep("group1", length(group1_samples)), rep("group2", length(group2_samples))))
    
    dds <- DESeqDataSetFromMatrix(countData = exp_data_filtered,
                                  colData = data.frame(condition = condition),
                                  design = ~ condition)
    
    dds <- DESeq(dds)
    res <- results(dds)
    
    log2_fold_changes <- res$log2FoldChange
    p_values <- res$pvalue
    
    if (adjust_pvalue) {
      adj_pvalues <- res$padj  # DESeq2 already provides adjusted p-values
    } else {
      adj_pvalues <- res$pvalue
    }
    
    deg_results <- data.frame(
      Symbols = rownames(res),
      baseMean = res$baseMean,
      Group1Mean = rowMeans(exp_data_filtered[, group1_samples], na.rm = TRUE),
      Group2Mean = rowMeans(exp_data_filtered[, group2_samples], na.rm = TRUE),
      Log2FoldChange = log2_fold_changes,
      lfcSE = res$lfcSE,
      stat = res$stat,
      PValue = res$pvalue,
      AdjPValue = adj_pvalues,
      stringsAsFactors = FALSE
    )
    
  } else {
    # Filter expression data to include only the relevant samples
    group1_data <- exp_data[, colnames(exp_data) %in% group1_samples, drop = FALSE]
    group2_data <- exp_data[, colnames(exp_data) %in% group2_samples, drop = FALSE]
    
    # Add pseudocount to avoid log transformation issues
    pseudocount <- 1
    group1_data <- group1_data + pseudocount
    group2_data <- group2_data + pseudocount
    
    # Ensure the data is numeric
    group1_data <- data.frame(lapply(group1_data, as.numeric))
    group2_data <- data.frame(lapply(group2_data, as.numeric))
    
    # Calculate log2 fold changes
    log2_fold_changes <- log2(rowMeans(group2_data, na.rm = TRUE)) - log2(rowMeans(group1_data, na.rm = TRUE))
    
    # Calculate p-values for each gene
    p_values <- vector("numeric", length = nrow(exp_data))
    for (i in 1:nrow(exp_data)) {
      gene_data1 <- as.numeric(group1_data[i, , drop = TRUE])
      gene_data2 <- as.numeric(group2_data[i, , drop = TRUE])
      
      if (test_type == "paired") {
        p_values[i] <- t.test(gene_data1, gene_data2, paired = TRUE)$p.value
      } else if (test_type == "unpaired") {
        p_values[i] <- t.test(gene_data1, gene_data2, paired = FALSE)$p.value
      } else if (test_type == "mann_whitney") {
        p_values[i] <- wilcox.test(gene_data1, gene_data2, exact = FALSE)$p.value
      } else if (test_type == "wilcoxon") {
        p_values[i] <- wilcox.test(gene_data1, gene_data2, paired = TRUE, exact = FALSE)$p.value
      }
    }
    
    # Adjust p-values if requested
    if (adjust_pvalue) {
      adj_pvalues <- p.adjust(p_values, method = "BH")
    } else {
      adj_pvalues <- p_values
    }
    
    deg_results <- data.frame(
      Symbols = rownames(exp_data),
      Group1Mean = rowMeans(group1_data, na.rm = TRUE),
      Group2Mean = rowMeans(group2_data, na.rm = TRUE),
      Log2FoldChange = log2_fold_changes,
      PValue = p_values,
      AdjPValue = adj_pvalues,
      stringsAsFactors = FALSE
    )
  }
  
  # Determine significance
  deg_results$Significant <- (deg_results$AdjPValue <= pvalue_cutoff) & (abs(deg_results$Log2FoldChange) >= log2(fc_cutoff))
  
  # Categorize regulation status
  deg_results$Regulation <- "Not Significant"
  deg_results$Regulation[deg_results$Significant & deg_results$Log2FoldChange >= log2(fc_cutoff)] <- "Upregulated"
  deg_results$Regulation[deg_results$Significant & deg_results$Log2FoldChange <= -log2(fc_cutoff)] <- "Downregulated"
  
  return(deg_results)
}


# Volcano plot generation function
generate_volcano_plot <- function(results, group1_name, group2_name) {
  total_genes <- nrow(results)
  up_genes <- sum(results$Regulation == "Upregulated", na.rm = TRUE)
  down_genes <- sum(results$Regulation == "Downregulated", na.rm = TRUE)
  
  p <- ggplot(results, aes(x = Log2FoldChange, y = -log10(PValue), color = Regulation,
                           text = paste("Gene:", Symbols,
                                        "<br>", group1_name, "Mean:", Group1Mean,
                                        "<br>", group2_name, "Mean:", Group2Mean,
                                        "<br>Fold Change:", Log2FoldChange,
                                        "<br>P-value:", PValue,
                                        "<br>Adj P-value:", AdjPValue))) +
    geom_point(alpha = 0.8) +
    labs(x = "Log2 Fold Change", y = "-log10(P-value)", title = "Volcano Plot") +
    scale_color_manual(values = c("Upregulated" = "red",
                                  "Downregulated" = "blue",
                                  "Not Significant" = "grey")) +
    theme_minimal() +
    coord_fixed()  
  
  ggplotly(p, tooltip = "text")
}




# MA plot generation function
generate_ma_plot <- function(exp_data, group1_samples, group2_samples, results, group1_name, group2_name) {
  # Add pseudocount to avoid log transformation issues
  pseudocount <- 1
  group1_data <- exp_data[, colnames(exp_data) %in% group1_samples, drop = FALSE] + pseudocount
  group2_data <- exp_data[, colnames(exp_data) %in% group2_samples, drop = FALSE] + pseudocount
  
  A <- (log2(rowMeans(group1_data)) + log2(rowMeans(group2_data))) / 2
  M <- log2(rowMeans(group2_data)) - log2(rowMeans(group1_data))
  
  ma_data <- data.frame(A = A, M = M,
                        Regulation = results$Regulation,
                        Symbols = results$Symbols,
                        Group1Mean = rowMeans(group1_data),
                        Group2Mean = rowMeans(group2_data),
                        Log2FoldChange = M,
                        PValue = results$PValue,
                        AdjPValue = results$AdjPValue)
  
  total_genes <- nrow(ma_data)
  up_genes <- sum(ma_data$Regulation == "Upregulated", na.rm = TRUE)
  down_genes <- sum(ma_data$Regulation == "Downregulated", na.rm = TRUE)
  
  print(paste("Total genes in MA plot:", total_genes))
  print(paste("Upregulated genes:", up_genes))
  print(paste("Downregulated genes:", down_genes))
  
  p <- ggplot(ma_data, aes(x = A, y = M, color = Regulation,
                           text = paste("Gene:", Symbols,
                                        "<br>", group1_name, "Mean:", Group1Mean,
                                        "<br>", group2_name, "Mean:", Group2Mean,
                                        "<br>Fold Change:", Log2FoldChange,
                                        "<br>P-value:", PValue,
                                        "<br>Adj P-value:", AdjPValue))) +
    geom_point(alpha = 0.5) +
    labs(x = "Average Log Intensity (A)", y = "Log Fold Change (M)", title = "MA Plot") +
    scale_color_manual(values = c("Upregulated" = "red",
                                  "Downregulated" = "blue",
                                  "Not Significant" = "grey")) +
    theme_minimal() +
    coord_fixed()
  
  ggplotly(p, tooltip = "text")
}




PairwiseComparisonTabServer <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Access the processed data from the dataset reactive
    expression_data <- reactive(dataset$filtered_data())
    groups_data <- reactive(dataset$groups_data())
    selected_groups <- reactive(dataset$selected_groups())
    
    # Reactive value to store the results
    deg_results <- reactiveVal()
    
    # Disable download buttons initially
    shinyjs::disable("download_all")
    shinyjs::disable("download_filtered")
    
    observe({
      req(selected_groups())
      grp_names <- selected_groups()
      if (is.null(grp_names)) return()
      updateSelectInput(session, "select_group1", choices = grp_names, selected = grp_names[1])
      updateSelectInput(session, "select_group2", choices = grp_names, selected = grp_names[length(grp_names)])
    })
    
    observeEvent(input$run_analysis, {
      req(expression_data(), groups_data())
      exp_data <- expression_data()
      grp_data <- groups_data()
      
      showModal(modalDialog("Running analysis, please wait...", footer = NULL))
      
      # Map group selections to sample names
      group1_samples <- grp_data$Sample[grp_data$Group == input$select_group1]
      group2_samples <- grp_data$Sample[grp_data$Group == input$select_group2]
      
      # Ensure the sample names match the column names in the expression data
      group1_samples <- intersect(group1_samples, colnames(exp_data))
      group2_samples <- intersect(group2_samples, colnames(exp_data))
      
      # Perform the DEG analysis
      res <- perform_DEG_analysis(exp_data, group1_samples, group2_samples,
                                  input$test_type, input$fc_cutoff, input$pvalue_cutoff, input$adjust_pvalue)
      
      # Store the results in a reactive value
      deg_results(res)
      
      # Enable download buttons after analysis
      shinyjs::enable("download_all")
      shinyjs::enable("download_filtered")
      
      # Separate results into positive and negative DEG tables
      pos_results <- res[res$Regulation == "Upregulated", ]
      neg_results <- res[res$Regulation == "Downregulated", ]
      
      # Calculate and display gene counts
      total_genes <- nrow(res)
      up_genes <- sum(res$Regulation == "Upregulated", na.rm = TRUE)
      down_genes <- sum(res$Regulation == "Downregulated", na.rm = TRUE)
      
      # Update gene counts
      output$volcano_gene_counts <- renderText({
        paste("Total genes:", total_genes,
              "| Upregulated genes:", up_genes,
              "| Downregulated genes:", down_genes)
      })
      
      output$ma_gene_counts <- renderText({
        paste("Total genes:", total_genes,
              "| Upregulated genes:", up_genes,
              "| Downregulated genes:", down_genes)
      })
      
      # Update outputs: Volcano Plot, MA Plot, DEG Tables
      output$volcano_plot <- renderPlotly({ 
        on.exit(removeModal())
        generate_volcano_plot(res, input$select_group1, input$select_group2) })
      
      output$ma_plot <- renderPlotly({ 
        on.exit(removeModal())
        generate_ma_plot(exp_data, group1_samples, group2_samples, res, input$select_group1, input$select_group2) })
      
      output$pos_deg_table <- renderDT({ datatable(pos_results, options = list(pageLength = 10), rownames = FALSE) })
      
      output$neg_deg_table <- renderDT({ datatable(neg_results, options = list(pageLength = 10), rownames = FALSE) })
      
    })
    
    # Download handlers for DEG results
    output$download_all <- downloadHandler(
      filename = function() { "full_deg_output.csv" },
      content = function(file) {
        req(deg_results())
        write.csv(deg_results(), file, row.names = FALSE)
      }
    )
    output$download_filtered <- downloadHandler(
      filename = function() { "filtered_deg_output.csv" },
      content = function(file) {
        req(deg_results())
        write.csv(subset(deg_results(), Regulation != "Not Significant"), file, row.names = FALSE)
      }
    )
  })
}
