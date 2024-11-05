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
                             plotlyOutput(ns("volcano_plot"))),
                    tabPanel("MA Plot",
                             textOutput(ns("ma_gene_counts")),
                             plotlyOutput(ns("ma_plot"))),
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
      adj_pvalues <- p.adjust(p_values, method = "BH")
    } else {
      adj_pvalues <- p_values
    }

    results <- data.frame(
      Symbols = exp_data$Symbol,
      Gene_Symbol = exp_data$Gene_Symbol,
      baseMean= res$baseMean,
      Group1Mean = rowMeans(exp_data_filtered[, group1_samples], na.rm = TRUE),
      Group2Mean = rowMeans(exp_data_filtered[, group2_samples], na.rm = TRUE),
      Log2FoldChange = log2_fold_changes,
      lfcSE = res$lfcSE,
      stat = res$stat,
      PValue = res$pvalue,
      padj= res$padj,
      Significant = (abs(log2_fold_changes) >= log2(fc_cutoff)) & (adj_pvalues <= pvalue_cutoff),
      stringsAsFactors = FALSE
    )

    return (results)

  }
  else {
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
    gene_data1 <- group1_data[i, , drop = FALSE]  # Ensure it remains as a dataframe
    gene_data2 <- group2_data[i, , drop = FALSE]  # Ensure it remains as a dataframe

    if (test_type == "paired") {
      p_values[i] <- t.test(as.numeric(gene_data1), as.numeric(gene_data2), paired = TRUE)$p.value
    } else if (test_type == "unpaired") {
      p_values[i] <- t.test(as.numeric(gene_data1), as.numeric(gene_data2), paired = FALSE)$p.value
    } else if (test_type == "mann_whitney") {
      p_values[i] <- wilcox.test(as.numeric(gene_data1), as.numeric(gene_data2), exact = FALSE)$p.value
    } else if (test_type == "wilcoxon") {
      p_values[i] <- wilcox.test(as.numeric(gene_data1), as.numeric(gene_data2), paired = TRUE, exact = FALSE)$p.value
    }
  }

  results <- data.frame(
    Symbols = exp_data$Symbol,
    Genes = exp_data$Gene_Symbol,
    Log2FoldChange = log2_fold_changes,
    PValue = p_values,
    Group1Mean = rowMeans(exp_data[, group1_samples], na.rm = TRUE),
    Group2Mean = rowMeans(exp_data[, group2_samples], na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  # Adjust p-values if requested
  if (adjust_pvalue) {
    results$AdjPValue <- p.adjust(results$PValue, method = "BH")
  } else {
    results$AdjPValue <- results$PValue
  }

  # Mark significant genes
  results$Significant <- (abs(results$Log2FoldChange) >= log2(fc_cutoff)) & (results$AdjPValue <= pvalue_cutoff)

  return(results)
  }
}

# Volcano plot generation function
generate_volcano_plot <- function(results, group1_name, group2_name) {
  total_genes <- nrow(results)
  red_genes <- sum(results$Significant, na.rm = TRUE)
  
  p <- ggplot(results, aes(x = Log2FoldChange, y = -log10(PValue), color = Significant,
                           text = paste("Gene:", Symbols,
                                        "<br>", group1_name, "Mean:", Group1Mean,
                                        "<br>", group2_name, "Mean:", Group2Mean,
                                        "<br>Fold Change:", Log2FoldChange,
                                        "<br>P-value:", PValue))) +
    geom_point(alpha = 0.8) +
    labs(x = "Log2 Fold Change", y = "-log10(P-value)", title = "Volcano Plot") +
    scale_color_manual(values = c("grey", "red"), na.value = "grey") +
    theme_minimal()
  
  print(paste("Total genes in volcano plot:", total_genes))
  print(paste("Red genes in volcano plot:", red_genes))
  
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
  
  ma_data <- data.frame(A = A, M = M, Significant = results$Significant,
                        Symbols = results$Symbols,
                        Group1Mean = rowMeans(group1_data),
                        Group2Mean = rowMeans(group2_data),
                        Log2FoldChange = M,
                        PValue = results$PValue)
  
  total_genes <- nrow(ma_data)
  red_genes <- sum(ma_data$Significant, na.rm = TRUE)
  
  print(paste("Total genes in MA plot:", total_genes))
  print(paste("Red genes in MA plot:", red_genes))
  
  p <- ggplot(ma_data, aes(x = A, y = M, color = Significant,
                           text = paste("Gene:", Symbols,
                                        "<br>", group1_name, "Mean:", Group1Mean,
                                        "<br>", group2_name, "Mean:", Group2Mean,
                                        "<br>Fold Change:", Log2FoldChange))) +
    geom_point(alpha = 0.5) +
    labs(x = "Average Log Intensity (A)", y = "Log Fold Change (M)", title = "MA Plot") +
    scale_color_manual(values = c("grey", "red"), na.value = "grey") +
    theme_minimal()
  
  ggplotly(p, tooltip = "text")
}



PairwiseComparisonTabServer <- function(id, dataset) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Access the processed data from the dataset reactive
    expression_data <- reactive(dataset$filtered_data())
    groups_data <- reactive(dataset$groups_data())
    selected_groups <- reactive(dataset$selected_groups())
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
      
      # Map group selections to sample names
      group1_samples <- grp_data$Sample[grp_data$Group == input$select_group1]
      group2_samples <- grp_data$Sample[grp_data$Group == input$select_group2]
      
      # Debug statements
      print(paste("Group 1 Samples:", toString(group1_samples)))
      print(paste("Group 2 Samples:", toString(group2_samples)))
      print(paste("Expression Data Columns:", toString(colnames(exp_data))))
      
      # Ensure the sample names match the column names in the expression data
      group1_samples <- intersect(group1_samples, colnames(exp_data))
      group2_samples <- intersect(group2_samples, colnames(exp_data))
      
      # Debug statements after intersection
      print(paste("Filtered Group 1 Samples:", toString(group1_samples)))
      print(paste("Filtered Group 2 Samples:", toString(group2_samples)))
      
      # Perform the DEG analysis
      results <- perform_DEG_analysis(exp_data, group1_samples, group2_samples,
                                      input$test_type, input$fc_cutoff, input$pvalue_cutoff, input$adjust_pvalue)
      
      # Get significant genes
      significant_genes <- results$Symbols[results$Significant]
      
      # Debug statement for significant genes
      print(paste("Significant Genes:", toString(significant_genes)))
      
      # Separate results into positive and negative DEG tables
      pos_results <- results[results$Log2FoldChange >= 0, ]
      neg_results <- results[results$Log2FoldChange < 0, ]
      
      # Calculate and display gene counts
      total_genes <- nrow(results)
      red_genes <- sum(results$Significant, na.rm = TRUE)
      
      output$volcano_gene_counts <- renderText({
        paste("Total genes:", total_genes, "Red genes:", red_genes)
      })
      
      output$ma_gene_counts <- renderText({
        paste("Total genes:", total_genes, "Red genes:", red_genes)
      })
      
      # Update outputs: Volcano Plot, MA Plot, DEG Tables
      output$volcano_plot <- renderPlotly({ generate_volcano_plot(results, input$select_group1, input$select_group2) })
      output$ma_plot <- renderPlotly({ generate_ma_plot(exp_data, group1_samples, group2_samples, results, input$select_group1, input$select_group2) })
      output$pos_deg_table <- renderDT({ datatable(pos_results, options = list(pageLength = 10), rownames = FALSE) })
      output$neg_deg_table <- renderDT({ datatable(neg_results, options = list(pageLength = 10), rownames = FALSE) })
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
